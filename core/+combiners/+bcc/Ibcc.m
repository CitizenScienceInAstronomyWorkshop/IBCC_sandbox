classdef  Ibcc < combiners.AbstractCombiner
 
    properties
        label
        
        %Experiment settings ----------------------------------------------
        
        settings = []; %a BccSettings object
                
        nScores       
                
        nKnown = 10;
        
        mapFunction = [];
        
        trustedAlpha = [];
                
        trainVoteCounts;
        trainIdxs;
        testIdxs;
        
        targetWeights = [];
        
        validDecs = [];
        noMissingScores = true; %if true, we assume that once the scoreset 
        %has been constructed, there are no "noScore" values so don't look for them.
        
        %Method Parameters ------------------------------------------------
        fixedAlpha = true;
        
        %Model Parameters -------------------------------------------------
        Alpha0 = [];
        Nu0;
        
        Pi; %final confusion matrix 
        lnPi;
        lnK;
        Alpha; %final posterior alphas to be returned as agent ratings
        C; %the last set of results seen
        Nu = [];
        
        Tmat; %targets in matrix form, columns with all-0s for unknown labels
        Cmat;
        CmatByScore;                
        Cl = {}; %binary vectors for each possible score indicating where they occurred in the scoreset
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
        function obj = Ibcc(bccSettings, nAgents, K, targets, ~, nClasses, nScores)
            obj@combiners.AbstractCombiner(nAgents, K, targets, nClasses);
            obj.settings = bccSettings;
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
                        
        function T = priorT(obj, P, nDataPoints)
            U = rand(1, nDataPoints);
            T = U > P(2);
        end
        
        %Find the p(Ti | Ci, P, Conf, Alpha) =  p(Ti | Ci, P, Conf)
        function postTi = posteriorTi(obj, C, Conf, P, nSamples, members)
            
            pCiTi = zeros(1, nSamples, obj.nClasses);
            for j=1:obj.nClasses % ti == j
                %rows - true classes
                %columns - classifier outputs
                %3rd dimension - classifiers
                
                indConfj = {};
                if ~obj.noMissingScores
                    validDecs = C{3}(C{3}~=obj.noScore);
%                 zeroMapDecs = zeros(1, length(C{3})) + obj.zeroMapScore;
%                 zeroMapDecs(validDecs~=0) = validDecs(validDecs~=0);
%                 validDecs = zeroMapDecs;
                    validA = C{1}(C{3}~=obj.noScore);
                    validN = C{2}(C{3}~=obj.noScore);
                else
                    validDecs = C{3};
                    validA = C{1};
                    validN = C{2};
                end
                indxs = sub2ind([size(Conf,2), size(Conf,3)], validDecs, validA);
                logPiValues = log(Conf(j,indxs));
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
            
            if ~obj.noMissingScores
                try
                    %why do we use the zero score map like this? If we are
                    %instantiating all the valid decisions that are given a
                    %score of zero, why not just put them in the scoreset C to
                    %start with? Adding obj.zeroScoreMap to all scores at
                    %each iteration is quicker than finding the specific indices 
                    %to correct, but it would be quicker to do this as a one-off addition.
                    nonZeroScores = C{3}~=obj.noScore;
                    nValid = sum(nonZeroScores);
%                      if isempty(obj.validDecs)
                        obj.validDecs = zeros(nValid, 1) + obj.settings.zeroScoreMap;
                        obj.validDecs(C{3}(nonZeroScores)~=0) = C{3}(C{3}~=0 & nonZeroScores);
%                      end
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
            else
                nValid = length(C{2});
                validN = C{2};
                realIdxs = unique(validN);
                realIdxMat = sparse(realIdxs, 1, ones(1,length(realIdxs)), nAssets, 1);
                indx = sub2ind([obj.nScores obj.nAgents], C{3}, C{1});                 
            end            
   
            lnpCT = sparse([], [], [], nAssets, obj.nClasses, length(C{1}));
            pT = sparse([], [], [], nAssets, obj.nClasses, length(C{1})*obj.nClasses);

            lnPiIndx = lnPi(:,indx);
            
            useLike = obj.settings.useLikelihood;
            
            parfor j=1:obj.nClasses
                lnpCT(:, j) = sparse(validN', ones(1,nValid), lnPiIndx(j,:), nAssets, 1);
                
                if ~useLike
                    lnpCT(:, j) = lnpCT(:, j) + realIdxMat.*lnKappa(j);
                end
            end

            expA = sparse(size(lnpCT,1), obj.nClasses);
            rescale = -min(lnpCT,[],2);            
            parfor j=1:obj.nClasses    

                expAj = exp(lnpCT(:,j)+ realIdxMat.*rescale);
                %stop any extreme values from causing Na
                expAj(expAj==Inf) = 1;
                expA(:,j) = expAj;
            end
            expB = sum(expA,2);

            %using loops so we avoid having to repmat with large matrices
            parfor j=1:obj.nClasses
                pT(:,j) = expA(:,j)./expB;
            end
            
            if length(obj.targets)<1
                ET = pT;
                return
            end
            
            ET = obj.Tmat;
            ET(:, obj.testIdxs) = pT(obj.testIdxs, :)';
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

                mergedCount = sum(Count(j,:,rareAgents),3);
                                scale = (agentCounts(rareAgents)./mergeLevel) .^ 2;
                Count(j, :, rareAgents) = Count(j,:,rareAgents) .* ...
                    repmat(scale,1,obj.nScores);
                
            end
        end
        
        function Count = voteCounts(obj, C, T, agentIdx)
            testVoteCounts = zeros(obj.nClasses, obj.nScores, obj.nAgents);    
            Count = zeros(obj.nClasses, obj.nScores, obj.nAgents);    
                        
            for j=1:obj.nClasses
                try
                    Tj = T(j,:)';
                    zeroScoreMap = obj.settings.zeroScoreMap;
                    CmatSlices = obj.CmatByScore;
                    for l=1:obj.nScores
                        if l==zeroScoreMap
                            continue;
                        end
                        
                        Count(j, l, :) = CmatSlices{l} * Tj;
                    end
                catch me
                    display(['voteCounts fail in IBCC ' me.message]);
                end
            end
%              Count = obj.mergeLowCounts(Count);
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
            else
                mappedC = mappedC - obj.settings.minScore + 1;
            end
            baseScoreSet{3} = mappedC;
            mappedC = baseScoreSet;           
            
            obj.CmatByScore = {};
            for l=1:obj.nScores
                obj.CmatByScore{l} = sparse(mappedC{1}, mappedC{2}, double(mappedC{3}==l), obj.nAgents, length(obj.targets));
            end
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

