classdef  WeightedSum < combiners.weighted.WeightedCombiner
    %PRECWEIGHTEDAVERAGING Averages the posteriors using weights derived
    %from precision of logodds. *This doesn't really make sense as there is
    %no reason to use precision of logodds to weight posteriors. Posteriors
    %are not Gaussians, but should be betas, so cannot find maximum 
    %likelihood by precision-weighted addition.*
    
    properties
        label = 'Weighted Sum';
        
        %combinedPost = [];
        %combinedLogodds = [];
        
        discreetInputs = false;

        onTheFly = true      
        roundCombination = false;
    end
    
    properties (Constant)
        shortLabel = 'weightedSum';
        shortLabelRounded = 'roundedWeightedSum';
    end
    
    methods
        function obj = WeightedSum(nAgents, K, targets, agents, windowSize)
            obj@combiners.weighted.WeightedCombiner(nAgents, K, targets, agents);
            if nargin > 4
                obj.windowSize = windowSize;
            end
        end
        
        function scores = getScores(obj, post)
            if obj.discreetInputs
                scores = round(post);
            else
                scores = post;
            end
            
            scores = (scores- obj.minScore) ./ (obj.maxScore-obj.minScore);
        end
        
        function W = updateAndCombine(obj, n, scores, post, T)
            if ~obj.sequential
                windowEnd = n+1;
            else
                windowEnd = n;
            end
            W = obj.weightCalculator.calculateWeights(scores, post, T, obj.combinedPost, windowEnd);
            %obj.combinationStep(post, invW, n);  %don't do this if we already know the true target label
            if obj.sequential
                obj.combinedPost = [obj.combinedPost T(1,n)-1];
            else
                obj.combinedPost = [];
            end
        end
        
        function combinedPost = combineDataPoints(obj, scores, W)
            
            totalWeights = sum(W, 1); 
            if size(W,2)==1
                W = repmat(W, 1, size(scores,2));
            end
                
            weightedVotes = sum(W.*scores, 1);
            
            %if some base classifiers do not give a score they are ignored.
            combinedPost = weightedVotes ./ totalWeights;
        end
        
        function combinedPost = combineScoreset(obj, post, T)
                        
            beta = 0.01;
            windowSize = 100;
            
            nSamples = max(post{2});
            
            [validAgents, ~, idxInValidAgents] = unique(post{1});
            nValidAgents = length(validAgents);
            
            [validSamples, ~, idxInValidSamples] = unique(post{2});
            nValidSamples = length(validSamples);            
            
            obj.createWeightCalculator([], nSamples);
            
            factors = zeros(nValidAgents, 1);
            prevFactors = zeros(nValidAgents, 1);
            
            nContributions = zeros(nValidAgents, 1);
%             contributions = sparse(nValidAgents, nValidSamples);           
            sContribs = cell(1,nValidSamples);
            aContribs = cell(1,nValidAgents);
            
            errors = cell(1,nValidAgents);%sparse(nValidAgents, nValidSamples);            
            
            scoreMatrix = sparse(idxInValidAgents, idxInValidSamples, post{3}, nValidAgents, nValidSamples);
            
            %go through each data point in turn and combine the answers. Then 
            %see if we have a target label. If so, update the weights.
            
            combinedPost = zeros(1, nValidSamples);
            
            for n=1:length(post{2})
                
                sample = idxInValidSamples(n);
                realSample = post{2}(n);
                
                currentAgent = idxInValidAgents(n);
                                           
%                 scoreMatrix(currentAgent, sample) = post{3}(n);
                
                %contribution toward combined result of each individual, i.e. how much they agree with the combined vote 
%                 contributions(currentAgent, sample) = n;
                sContribs{sample} = [sContribs{sample} currentAgent];
                aContribs{currentAgent} = [aContribs{currentAgent} sample];
                contributors = sContribs{sample};%find(contributions(:,sample)~=0);
                
                %combine scores for this data point with non-updated
                %weights but using the new information from currentAgent

                nMissingData = windowSize-nContributions(contributors);
                nMissingData(nMissingData<0) = 0;
                oldFactors = factors(contributors);
                powerTerm = oldFactors+(nMissingData*0.5);
                normFactor = min(powerTerm);
                Wold = beta .^ (powerTerm-normFactor);
                scores = obj.getScores(scoreMatrix(contributors,sample));
                combinedSample = obj.combineDataPoints(scores, Wold);                
                
                if isnan(factors(currentAgent))
                    display('WS error');
                end
                
                if T(realSample) == 0
                    continue;
                end

                normTarget = T(realSample) - 1;
                if normTarget < 0
                    normTarget = combinedSample;
                end
                
                nContributions(currentAgent) = nContributions(currentAgent) + 1;
                
                normScores = (scoreMatrix(currentAgent,sample) - obj.minScore) ./ (obj.maxScore-obj.minScore);
                
                contribution = abs(normTarget-normScores) .* abs(normTarget-combinedSample);

%                 errors(currentAgent, sample) = contribution;       
                errors{currentAgent} = [errors{currentAgent} contribution];

                factors(currentAgent) = prevFactors(currentAgent) + contribution;

                if nContributions(currentAgent) > windowSize
                    expiredIdx = aContribs{currentAgent}(1); %find(contributions(currentAgent,:), 1);
                    expiredContrib = errors{currentAgent}(1);% expiredIdx);
                    factors(currentAgent) = factors(currentAgent) - expiredContrib;
                    errors{currentAgent}(1) = [];
%                     contributions(currentAgent, expiredIdx) = 0;
                    aContribs{currentAgent}(1) = [];
                    sIdx = sContribs{expiredIdx}==currentAgent;
                    sContribs{expiredIdx}(sIdx) = [];
                end

                normFactor = min(factors(factors~=0));
                if isempty(normFactor)
                    normFactor = 0;
                end
                Wnew = beta .^ (factors(contributors)-normFactor);
                scores = obj.getScores(scoreMatrix(contributors,sample));  
                combinedPost(sample) = obj.combineDataPoints(scores, Wnew);
                
                if combinedPost(sample)<0 || combinedPost(sample) > 1 || isnan(combinedPost(sample)) 
                    display(['WS has invalid value: ' num2str(combinedPost(sample))]);
                end
                
                prevFactors = factors;
                
                if mod(n,5000)==0
                    display(['weighted sum sample ' num2str(n) ' of ' num2str(length(post{2}))]);
                end
            end      
            
            combinedPost = sparse(1, validSamples, combinedPost, 1, nSamples);
        end
        
        function setSequential(obj, seq)
            obj.sequential = true;
            obj.useFinalWeights = seq;
        end
                
        function combinedPost = combineCluster(obj, post, T, members, k)                   
            
            if obj.scoreSet
                combinedPost = obj.combineScoreset(post, T);
                return;
            end
            
            scores = obj.getScores(post);            
            
            %this may need to know which agents are in the cluster!
            obj.createWeightCalculator(members, size(post, 2));
                        
            obj.combinedPost = [];
            
            nTargets = length(T);
            
            if nTargets==0
                nTargets = 1;
            end
            
            W = [];
            if obj.sequential
                startIdx = find(T~=-1) + 1;
                startIdx = [startIdx nTargets+1];
            else
                startIdx = nTargets;
            end
            
            nLearnt = 1;
                        
            Trep = repmat(T, obj.nAgents, 1);
            
            
            if obj.sequential
                n = 2;
                W = obj.updateAndCombine(1, scores, post, Trep);
            else
                W = obj.updateAndCombine(nTargets, scores, post, Trep);
                obj.combinationStep(post, W, 1, nTargets); 
                n = nTargets + 1; %skip the next loop
            end
            while n <= nTargets
                if n < startIdx(1)
                    obj.combinationStep(post, W(:,end), n, startIdx(1)-1); 
                    n = startIdx(1);
                elseif n==startIdx(1);
                    W = obj.updateAndCombine(n, scores, post, Trep);
                    n = n+1;
                    nLearnt = nLearnt+1;
                    startIdx = startIdx(2:end);
                end
            end
            
            
            if obj.useFinalWeights
                obj.combinedPost = [];
                obj.combinationStep(post, W(:,end), 1, nTargets);
            end

            combinedPost = obj.combinedPost;
        end
        
        function plotWeights(obj, W, k)
            figure;
            
            labels = {size(W, 1)};
            for a=1:size(W, 1)
                plot(W(a,:), 'LineWidth', mod(a,3)+1);
                hold all
                
                labels{a} = sprintf('weight for agent %d', a); 
            end
            title(sprintf('%s, cluster %d', obj.label, k));
            
            legend(labels);
            hold off
        end
        
        function createWeightCalculator(obj, members, nSamples)
            if nargin < 2
                members = [];
            end
            obj.weightCalculator = combiners.weighted.weights.CombinedErrorWeights(obj.nAgents, obj.windowSize, members, nSamples);
        end
        
        function baseData = correctBaseData(obj, baseOutputs)
            
            if iscell(baseOutputs) && (length(baseOutputs)==3 || length(baseOutputs)==4)
                baseData = baseOutputs;
                if length(baseOutputs)==4
                    obj.targets = baseOutputs{4};
                end
                obj.scoreSet = true;
                obj.nAgents = max(baseData{1});                
            elseif obj.useSparse
                baseData = cell(1, 3);
                baseData{3} = reshape(baseOutputs, numel(baseOutputs), 1);
                baseData{1} = reshape(...
                    (1:size(baseOutputs,1))' * ones(1, size(baseOutputs,2)), ...
                    numel(baseOutputs), 1);
                baseData{2} = reshape(...
                    ones(size(baseOutputs,1), 1) * (1:size(baseOutputs,2)), ...
                    numel(baseOutputs), 1);                
                obj.scoreSet = true;
                obj.nAgents = max(baseData{1});                
            else
                obj.nAgents = size(baseOutputs, 1);
                baseData = baseOutputs;
            end        
        end        
        
        function combinedPost_n = combinationStep(obj, post, W, n, nEnd)
                        
            if nargin < 5
                nEnd = n;
            end
            
            %use the weights to weight the posteriors                       
            if n > size(W, 2) || nEnd > n
                W = W(:, end);
                post = post(:, n:nEnd);
            else
                post = post(:, n);
                W = W(:, n);
            end
%             %normalise so that smallest weight will be one
%             nW = max(invW) ./ invW;
%             
%             totalWeights = sum(nW); 
%             
%             weightedVotes = nW' * post;
            
            %if some base classifiers do not give a score they are ignored.
            combinedPost_n = obj.combineDataPoints(post, W);%weightedVotes ./ totalWeights;
            
            obj.combinedPost = [obj.combinedPost combinedPost_n];
        end
    end
end

