classdef  Ibcc < combiners.AbstractCombiner
 
    properties
        
        settings = []; %a BccSettings object
        
        agents

        Nu0;
        
        label
        
        nScores

        dIdxs
        
        nKnown = 10;
        
        samplingType = 'Gibbs'; %or 'MH';      
                
        minConf = 10 ^ -6;
        
        %maximum number of times we will run the gibbs sampler before
        %giving up
        maxSampleAttempts = 20;
        
        confSave = '';
        
        datasetLabel = '';
        
        fixedAlpha = true;
        Alpha0 = [];
        mapFunction = [];
        
        trustedAlpha = [];
        
        scoreSet = true;
        
        Pi; %final confusion matrix 
        lnPi;
        lnK;
        Alpha; %final posterior alphas to be returned as agent ratings
        C; %the last set of results seen
        Nu = [];
        
        trainVoteCounts;
        trainIdxs;
        testIdxs;
        Tmat; %targets in matrix form, columns with all-0s for unknown labels
        Cmat;

        targetWeights = [];
                
        Cl = {}; %binary vectors for each possible score indicating where they occurred in the scoreset
        validDecs = [];
    end
    
    properties (Access=protected)
        trustFinalAgent = [];        
    end
    
    methods (Static)        
        function f = dirpdf(X, A)
            B = prod(gamma(A), 2) ./ gamma(sum(A, 2));
            
            f = 1./B .* prod(X.^(A-1), 2);
        end
        
        function f=logDirPdf(logX, A)
            invB = log(gamma(sum(A, 2))) - sum(log(gamma(A)), 2);
            
            f = invB + sum(logX.*(A-1), 2);           
        end
    end
    
    methods       
        function obj = Ibcc(bccSettings, nAgents, K, targets, agents, nClasses, nScores)
            obj@combiners.AbstractCombiner(nAgents, K, targets, nClasses);
            obj.settings = bccSettings;
            obj.agents = agents;
            obj.normalised = true;
                        
            obj.label = 'I.B.C.C.';
            
            %obj.noScore = 0;
            
            obj.debug = obj.settings.debug;

            if nargin > 4
                obj.nClasses = nClasses;
                obj.nScores = nClasses;
            else
                obj.nClasses = 2;
                obj.nScores = 2;
            end   
            
            if nargin > 5
                obj.nScores = nScores;
            end
            
            obj.Nu0 = bccSettings.nu{bccSettings.iNu};
            if length(obj.Nu0)<obj.nClasses
                obj.Nu0 = [obj.Nu0 ones(1,obj.nClasses-length(obj.Nu0)).*obj.Nu0(1)];
            end
            
            obj.setAlphaPrior(bccSettings.Alpha);
                                    
            obj.combinerInfo = ['nu ' mat2str(obj.Nu0)];            
            obj.nonDet = true; %non-deterministic
        end
        
        function setTrustedAgents(obj, trustFinalAgent)
            if nargin > 1
                obj.trustFinalAgent = trustFinalAgent;
            end
        end
        
        function setTargets(obj, targets)
            obj.targets = targets;
            obj.Tmat = zeros(obj.nClasses, length(targets));
            obj.trainIdxs = find(targets~=0);
            obj.testIdxs = find(targets==0);
            if isempty(obj.targetWeights)
                obj.Tmat(sub2ind([size(obj.Tmat,1) size(obj.Tmat,2)], ...
                    targets(obj.trainIdxs), obj.trainIdxs)) = 1;
            else
                 obj.Tmat(sub2ind([size(obj.Tmat,1) size(obj.Tmat,2)], ...
                    targets(obj.trainIdxs), obj.trainIdxs)) = obj.targetWeights(obj.trainIdxs);               
            end
        end        
                
        %use an explicit alpha prior instead of the lambda stuff
        function setAlphaPrior(obj, AlphaPrior, dupe, nDupes)
            obj.Alpha0 = AlphaPrior;
            
            if nargin>3 && dupe && size(AlphaPrior,3)<nDupes        
                if size(AlphaPrior,3)==1            
                    obj.Alpha0 = repmat(obj.Alpha0, 1, nDupes);
                    obj.Alpha0 = reshape(obj.Alpha0, size(obj.Alpha0,1), obj.nScores, nDupes);
                else
                    obj.Alpha0 = zeros(size(obj.Alpha0,1), obj.nScores, nDupes);
                    obj.Alpha0(:,:,1:size(AlphaPrior,3)) = AlphaPrior;
                    for a=size(AlphaPrior,3)+1:nDupes
                        obj.Alpha0(:,:,a) = ...
                            sum(AlphaPrior,3) ./ size(AlphaPrior,3);
                    end
                end    
                for trusted=obj.trustFinalAgent
                    trustedAgent = obj.nAgents-trusted+1;
                    obj.Alpha0(:,:,trustedAgent) = obj.trustedAlpha(:,:,trusted);
                end
            elseif nargin > 3 && size(AlphaPrior,3)>nDupes
                obj.Alpha0 = AlphaPrior(:,:,1:nDupes);
            end
        end
                
        function Alpha = priorAlpha(obj)
            
            if ~isempty(obj.Alpha0)
                Alpha = obj.Alpha0;
                return
            end
            
            if obj.fixedAlpha
                Alpha = 1 ./ obj.lambda;
            else
                %support for confusion matrix. Pick the mode of the priors.
                Alpha = ones(obj.nClasses, obj.nScores, obj.nAgents);
                for i=1:obj.nAgents
                    Alpha(:, :, i) = exprnd(1./obj.lambda(:,:,i), obj.nClasses, obj.nScores);
                end
            end
            obj.Alpha0 = Alpha;            
        end
        
        function P = priorP(obj)
            %class proportions. 
            P = gamrnd(obj.Nu0, ones(1, obj.nClasses));
            P = P ./ (sum(P)*ones(1, obj.nClasses));
        end        
        
        function ConfSample = sampleConf(obj, Alpha, C)
            agents = unique(C{1});

            ConfSample = zeros(obj.nClasses, obj.nScores, obj.nAgents);
            
            if obj.nScores > 2 || obj.nClasses > 2
                for j=1:obj.nClasses
                    K = agents;
                    while ~isempty(K)
                        for c=1:obj.nScores
                            ConfSample(j, c, K) = betarnd(reshape(Alpha(j,c,K), 1, length(K)), ...
                                reshape(sum(Alpha(j,:,K),2)-Alpha(j,c,K), 1, length(K)) );

                        end
                        [J, K] = find(ConfSample(j, :, agents)<obj.minConf, 1);    
                        ConfSample(j, K) = obj.minConf;
                        K = [];
                    end
                end
            else
                for j=1:obj.nClasses
                    K = 1:obj.nAgents;
                    while ~isempty(K)
                        ConfSample(j, 1, K) = betarnd(Alpha(j, 1, K), Alpha(j, 2, K));
                        ConfSample(j, 2, K) = 1 - ConfSample(j, 1, K);
                        [J, K] = find(ConfSample(j, :, :)<obj.minConf, 1);
                        K = ceil(K ./ obj.nClasses);
                    end
                end             
            end
        end  
                
        function T = priorT(obj, P, nDataPoints)
            U = rand(1, nDataPoints);
            T = U > P(2);
        end
        
        %THIS ONLY WORKS WITH TWO CLASSES
        function [T, postTi] = sampleT(obj, C, Conf, P,  members)     
            nSamples = length(obj.targets);
            
            %display('breakage');
            %C = C(:, length(Tknown)+1:nSamples);
            
            if nargin > 5 && ~isempty(members)
                display('breakage breakage breakage');
                %members = members(:, length(Tknown)+1:nSamples);
                postTi = obj.posteriorTi(C, Conf, P, nSamples, members);
            else
                postTi = obj.posteriorTi(C, Conf, P, nSamples);
            end
            
            T = zeros(1, length(postTi));
            if obj.scoreSet
                validT = unique(C{2});
            else
                validT = 1:length(postTi);
            end
            randoms = rand(1, length(validT));
            
            T(validT) = (randoms<postTi(validT)) + 1;
            
            %display('breaking tknowns again');
            T(obj.targets~=0) = obj.targets(obj.targets~=0);
            
            %T = [obj.targets Tsample];
        end
        
        %Find the p(Ti | Ci, P, Conf, Alpha) =  p(Ti | Ci, P, Conf)
        function postTi = posteriorTi(obj, C, Conf, P, nSamples, members)
            
            pCiTi = zeros(1, nSamples, obj.nClasses);
            for j=1:obj.nClasses % ti == j
                %rows - true classes
                %columns - classifier outputs
                %3rd dimension - classifiers
                
                indConfj = {};
                
                validDecs = C{3}(C{3}~=obj.noScore);
                zeroMapDecs = zeros(1, length(C{3})) + obj.zeroMapScore;
                zeroMapDecs(validDecs~=0) = validDecs(validDecs~=0);
                validDecs = zeroMapDecs;
                validA = C{1}(C{3}~=obj.noScore);
                validN = C{2}(C{3}~=obj.noScore);

                logPiValues = log(Conf(sub2ind(size(Conf), ...
                    j.*ones(length(validDecs), 1), validDecs, validA)));
                pCi_giv_Ti_eqj = sparse(validA, validN, logPiValues, obj.nAgents, nSamples);
                pCi_giv_Ti_eqj = exp(sum(pCi_giv_Ti_eqj, 1));
                
                %rows in pCi_given_Ti correspond to true labels. Columns
                %correspond to samples.
                if nargin > 5
                    %reject any non-members of the cluster by setting them
                    %to 1
                    pCi_giv_Ti_eqj(members==0) = 1;
                end
                
                pCiTi(1, :, j) = prod(pCi_giv_Ti_eqj, 1) .* P(j);
            end
            
            %the probabilities we calculate here also assume certain values
            %of hyperparameters - these conditions are left out of variable names for
            %simplicity. 
                       
            %return the probability that Ti==1, using the labelling scheme
            %[0 1].
            postTi = pCiTi(:, :, 2) ./ sum(pCiTi, 3);
            
            %multiply by p of the hyperparameters later or integrate over
            %them (not implemented yet). If we integrate over them to
            %produce P, conf and alpha we don't need to include them as a 
            %condition in the posterior.
        end
        
        function [ET, pT, lnpCT] = expectedT(obj,  C, lnKappa, lnPi, nAssets)     
            
            if nargin < 5
                nAssets = length(obj.targets);
            end
            
            try
                %why do we use the zero score map like this? If we are
                %instantiating all the valid decisions that are given a
                %score of zero, why not just put them in the scoreset C to
                %start with? Adding obj.zeroScoreMap to all scores at
                %each iteration is quicker than finding the specific indices 
                %to correct, but it would be quicker to do this as a one-off addition.
                nonZeroScores = C{3}~=obj.noScore;
                nValid = sum(nonZeroScores);
                if isempty(obj.validDecs)
                    obj.validDecs = zeros(nValid, 1) + obj.settings.zeroScoreMap;
                    obj.validDecs(C{3}(nonZeroScores)~=0) = C{3}(C{3}~=0 & nonZeroScores);
                end
                validN = C{2}(nonZeroScores);
                validAgents = C{1}(nonZeroScores);
                realIdxs = unique(validN);  
                
                realIdxMat = sparse(realIdxs, 1, ones(1,length(realIdxs)), nAssets, 1);
                
            catch me
                display(['index for E[T] failed: ' me.message]);
            end
            
            indx = sub2ind([obj.nScores obj.nAgents], ...
                            obj.validDecs, ...
                            validAgents);     
            lnpCT = sparse([], [], [], nAssets, obj.nClasses, length(C{1}));
            pT = sparse([], [], [], nAssets, obj.nClasses, length(C{1})*obj.nClasses);

            for j=1:obj.nClasses
                lnpCT(:, j) = sparse(validN', ones(1,nValid), lnPi(j, indx), nAssets, 1);
                
                if ~obj.settings.useLikelihood
                    lnpCT(:, j) = lnpCT(:, j) + realIdxMat.*lnKappa(j);
                end

%                 for i=1:nAssets
%                     for s=1:obj.nScores
%                         lnpCT(i,j) = sum(lnPi(j,s, obj.Cmat(:,i)==s));
%                     end
%                     lnpCT(i,j) = lnpCT(i,j) + lnKappa(j);
%                 end
            end

            expA = sparse(size(lnpCT,1), obj.nClasses);
            for j=1:obj.nClasses
                rescale = -min(lnpCT(:, :),[],2);

                expAj = exp(lnpCT(:,j)+ realIdxMat.*rescale);
                %stop any extreme values from causing Na
                expAj(expAj==Inf) = 1;
                expA(:,j) = expAj;
            end
            expB = sum(expA,2);

            %using loops so we avoid having to repmat with large matrices
            for j=1:obj.nClasses
                pT(:,j) = expA(:,j)./expB;
            end
            %data points that have no base classifier responses (stuck at
            %the priors).
%             unsetIdxs = obj.testIdxs(~ismember(obj.testIdxs,realIdxs));
%             nUnset = length(unsetIdxs);
%             for j=1:obj.nClasses
%                 unsetVal = repmat( exp(lnKappa(j)) ./ sum(exp(lnKappa)), nUnset, 1);
%                 pT(unsetIdxs, j) = unsetVal;
%             end
            
            if length(obj.targets)<1
                ET = pT;
                return
            end
            
            ET = obj.Tmat;
            ET(:, obj.testIdxs) = pT(obj.testIdxs, :)';
            
%             ET = pT;
%             ET(:, obj.targets~=0) = 0;
%             ET(sub2ind(size(ET), double(obj.targets(obj.targets~=0)), find(obj.targets~=0) )) = 1;
%             ET(1, obj.targets~=0) = 1-ET(2, obj.targets~=0);
        end     
        
        function Count = mergeLowCounts(obj,Count,mergeLevel)
            if nargin < 3
                mergeLevel = 0;
            end
            
            for j=1:obj.nClasses
                agentCounts = sum(Count(j,:,:),2);
                rareAgents = find(agentCounts>0 & agentCounts<mergeLevel);
                
                if isempty(rareAgents)
%                     display('no low counts to merge');
                    continue;
                end
                
%                 display(['**merging at level ' num2str(mergeLevel) ': ' num2str(length(rareAgents)) ' of total agents ' num2str(length(find(agentCounts>0)))]);
                
%                 maxCount = max(agentCounts);
                
                mergedCount = sum(Count(j,:,rareAgents),3);
                
%                 Count(j,:,rareAgents) = repmat(mergedCount, [1 1 length(rareAgents)]);
                scale = (agentCounts(rareAgents)./mergeLevel) .^ 2;
                Count(j, :, rareAgents) = Count(j,:,rareAgents) .* ...
                    repmat(scale,1,obj.nScores);
                
            end
            
            
%             agentCounts = full(sum(obj.Cmat~=obj.noScore,2));
%             rareAgents = find(agentCounts>0 & agentCounts<10);
%             if isempty(rareAgents)
%                 return;
%             end
%             rareId = rareAgents(1);
%             idxsToModify = ismember(C{1}, rareAgents);
%             
%             display(['IBCCVB: merging little-known agents: ' num2str(length(rareId)) ' affecting ' num2str(length(idxsToModify)) ' responses.']);
%             
%             C{1}(idxsToModify) = rareId;
%             respCounts = sparse(C{1}, C{2}, 1, obj.nAgents, length(obj.targets));
%             C3_avg = C{3} ./ respCounts(sub2ind(size(respCounts),C{1},C{2}));
%             %C3_avg = round(C3_avg);
%             obj.Cmat = sparse(C{1}, C{2}, C3_avg, obj.nAgents, length(obj.targets));
        end
        
        function Count = voteCounts(obj, C, T, agentIdx)
            testVoteCounts = zeros(obj.nClasses, obj.nScores, obj.nAgents);    
            Count = zeros(obj.nClasses, obj.nScores, obj.nAgents);    
                        
            if isempty(obj.trainVoteCounts)
                display(['no. agents: ' num2str(obj.nAgents) ]);
                
                trainScores = find(obj.targets(C{2})~=0);
                scores = C{3}(trainScores);
                
                display(['no. training points: ' num2str(length(trainScores)) ]);
                display(['highest score: ' num2str(max(scores)) ]);
                
                
                if max(scores)<obj.nScores
                    %we assume the zeros are meant to be the missing
                    %classifier output value
                    adjScores = zeros(length(scores), 1) + obj.settings.zeroScoreMap;
                    adjScores(scores~=0) = scores(scores~=0);
                    display('adjusting scores...');
                else
                    %we assume the zeros are intentionally skipping that
                    %classifier output
                    adjScores = scores;
                end
                
                display(['lowest score: ' num2str(min(adjScores)) ]);
                display(['highest agent idx: ' num2str(max(C{1}(trainScores))) ]);
                display(['number of targets: ' num2str(size(T,2)) ]);
                display(['highest target idx: ' num2str(max(C{2}(trainScores))) ]);
                                
                for j=1:obj.nClasses
                    %create the sparse with a small, fixed number of
                    %columns for speed, then transpose.
%                     Countj = sparse(C{1}(trainScores), adjScores, T(j,
%                     C{2}(trainScores)), obj.nAgents, obj.nScores);
%                     Count(j, :, :) = Countj';
                    for s=1:obj.nScores
                        Count(j, s, :) = (obj.Cmat==s) * T(j,:)';%sum(repmat(T(j,:), obj.nAgents, 1) .* obj.Cmat==s, 2);
                    end
                end
                obj.trainVoteCounts = Count;
            end
%             testScores = find(obj.targets(C{2})==0);            
%             C{1} = C{1}(testScores);
%             C{2} = C{2}(testScores);
%             C{3} = C{3}(testScores);
%             
%             if obj.scoreSet && isempty(C{1})
%                  return
%             end
%             
%             if nargin > 3
%                 C{3} = C{3}(C{1}==agentIdx);
%                 C{2} = C{2}(C{1}==agentIdx);
%                 C{1} = ones(length(C{2}),1);
%                 nAgents = 1;
%             else
%                 nAgents = obj.nAgents;
%             end
            
%             scores = C{3};
%             adjScores = zeros(length(C{3}), 1) + obj.zeroScoreMap;
%             adjScores(scores~=0) = scores(scores~=0);
            
            for j=1:obj.nClasses
                try
                    for l=1:obj.nScores
                        if l==obj.settings.zeroScoreMap
                            continue;
                        end
                        
%                         if l>length(obj.Cl) || isempty(obj.Cl{l})
%                             obj.Cl{l} = (C{3}==l)';
%                         end
                        Count(j, l, :) = (obj.Cmat==l) * T(j,:)';
%                         testVoteCounts(j,l,:) = sparse(C{1}, 1, T(j, C{2}) .* obj.Cl{l}, nAgents, 1);      
%                         Count(j, l, :) = testVoteCounts(j, l, :) + obj.trainVoteCounts(j, l, :);
                    end
                    if obj.settings.zeroScoreMap > 0
                        testVoteCounts(j, obj.settings.zeroScoreMap, :) = sum(T(j,obj.testIdxs)) - sum(testVoteCounts(j,:,:),2);                    
                        Count(j, obj.settings.zeroScoreMap, :) = testVoteCounts(j, obj.settings.zeroScoreMap, :) + obj.trainVoteCounts(j, obj.settings.settings.zeroScoreMap, :);
                    end
                    
%                     Count(j, :, :) = sparse(adjScores, C{1}, T(j, C{2}), obj.nScores, obj.nAgents);
%                     Count(j, :, :) = Count(j, :, :) + obj.trainVoteCounts(j, :, :);

%                     for i=1:length(C{1})
%                         if C{3}(i)==0 && C{3}(i)~=obj.noScore
%                             c = obj.zeroScoreMap;
%                         else
%                             c = C{3}(i);
%                         end
%                         Count(j, c, C{1}(i)) = Count(j, c, C{1}(i)) + T(j, C{2}(i));
%                     end
                catch me
                    display(['voteCounts fail in IBCC ' me.message]);
                end
            end
            
            Count = obj.mergeLowCounts(Count);
        end
        
        function mappedC = prepareC(obj, C)
            
            baseScoreSet = C;
            C = baseScoreSet{3};
            try
                C = round(C);
            catch me
                display('IBCC failed to round base classifier scores');
            end
            mappedC = C;
            
            if ~isempty(obj.settings.IbccMapFunction)
                mappedC = obj.settings.IbccMapFunction(C);
            elseif ~isempty(obj.settings.scoreMap)

                for s=1:size(obj.settings.scoreMap,2)
                    scoreIdxs = C==obj.settings.scoreMap(2,s);
                    mappedC(scoreIdxs) = obj.settings.scoreMap(1, s);
                end

                [minVal, minIdx] = min(obj.settings.scoreMap(2,:));
                [maxVal, maxIdx] = max(obj.settings.scoreMap(2,:));
                mappedC(C<minVal) = obj.settings.scoreMap(1,minIdx);
                mappedC(C>maxVal) = obj.settings.scoreMap(1,maxIdx);
            end
            baseScoreSet{3} = mappedC;
            mappedC = baseScoreSet;
        end

        function [combinedPost, agentRatings] = combineDecisions(obj, baseOutputs, clusters)

            nSamples = length(obj.targets);
            
            combinedPost = zeros(obj.K, nSamples);
            obj.Alpha = zeros(obj.nClasses, obj.nScores, obj.nAgents);            
            
            %get the relevant inputs
            if obj.K == 1
                combinedPost(1, :) = ...
                    obj.combineCluster( baseOutputs, ...
                        obj.targets);    
            else
                display('*** ignoring cluster of stupid agents: starting at k=2');
                for k=2:obj.K
                    %clusters change over time, so dlr needs to take different
                    %inputs over time - this is not possible.
                    %So we set values for base learners outside the cluster
                    %to 0.5 (regardless of label) as we can't ignore them.
                    
                    inputPosts = baseOutputs;
                    inputPosts(clusters~=k) = 0.5;
                                        
                    combinedPost(k, :) = ...
                        obj.combineCluster( inputPosts, ...
                            obj.targets, (clusters==k));
                end
            end
            agentRatings = obj.Alpha;
        end
        
        function baseData = correctBaseData(obj, baseOutputs)
            if iscell(baseOutputs) && (length(baseOutputs)==3 || length(baseOutputs)==4)
                baseData = baseOutputs;
                if length(baseOutputs)==4
%                     obj.targets = baseOutputs{4};
                    targets = baseOutputs{4};
                    labelIdxs = find(targets ~= 0);
                    agentIdxs = (obj.nAgents) * ones(length(labelIdxs),1);
                    decs = targets(labelIdxs);
                    
                    baseData{1} = [agentIdxs; baseData{1}];
                    baseData{2} = [labelIdxs'; baseData{2}];
                    baseData{3} = [decs'; baseData{3}];
                end
            else
                baseData = cell(1, 3);
                baseData{3} = reshape(baseOutputs, numel(baseOutputs), 1);
                baseData{1} = reshape(...
                    (1:size(baseOutputs,1))' * ones(1, size(baseOutputs,2)), ...
                    numel(baseOutputs), 1);
                baseData{2} = reshape(...
                    ones(size(baseOutputs,1), 1) * (1:size(baseOutputs,2)), ...
                    numel(baseOutputs), 1);                
            end 
            
            obj.scoreSet = true;
            obj.nAgents = max(baseData{1});
        end
    
        function savePi(obj, lnPi, C, nIt)
            
            realAgents = unique(C{1});
            Pi = exp(lnPi);
            Pi = Pi(:,:,realAgents);
                        
            save(sprintf('/homes/49/edwin/matlab/combination/data/galaxyZoo3/confusion_matrix/%dFold/%d.mat', obj.fold, nIt), 'Pi');
            
            if nIt==0
                agents = C{1};
                assets = C{2};
                decisions = C{3};
                save(sprintf('/homes/49/edwin/matlab/combination/data/galaxyZoo3/confusion_matrix/%dFold/agents.mat', obj.fold), 'agents');
                save(sprintf('/homes/49/edwin/matlab/combination/data/galaxyZoo3/confusion_matrix/%dFold/assets.mat', obj.fold), 'assets');
                save(sprintf('/homes/49/edwin/matlab/combination/data/galaxyZoo3/confusion_matrix/%dFold/decisions.mat', obj.fold), 'decisions');
            end
        end    
    end    
end

