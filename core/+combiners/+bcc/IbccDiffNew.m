classdef  IbccDiff < combiners.bcc.IbccVb
    %Use two conf matrix rows per class -> allows for different behaviour with
    %"difficult" items. Problem with this is we need to make assumptions apriori
    %about behaviour with different difficulty levels. Perhaps a better
    %alternative would be to perform unsupervised clustering of behavioural
    %patterns so that we could detect an array of latent sub-classes of items.
 
    properties 
	%blending parameter that controls the mixing between the data for a row in a confusion matrix 
	%and the pooled vector; it also controls the way that the pooled vector 
	%weights data from each row in the matrix
	B = 15;%10;%15;%value used in the emailed results %referred to as alpha in Tom Leonard 1977 renamed to avoid confusion.
	
	beta;
	gamma;
    end
 
    methods (Static)
        function sl = shortLabel
            sl = 'IbccDiff';
        end   
    end
    
    methods       
        function obj = IbccDiff(bccSettings, nAgents, K, targets, agents, nClasses, nScores)
            obj@combiners.bcc.IbccVb(bccSettings, nAgents, K, targets, agents, nClasses, nScores);
            obj.label = 'VB I.B.C.C. with pooling';
            
            obj.beta = 0.1;
            obj.gamma = repmat([1/3 1/3 1/3], [1, 1, obj.nAgents]);
            obj.settings.priorAdjLevel = 0; %ignore this bit of fiddling - this class is intended to avoid it completely
            
            obj.Nu0 = bccSettings.nu{bccSettings.iNu};
            if length(obj.Nu0)<obj.nClasses
                obj.Nu0 = [obj.Nu0 ones(1,obj.nClasses-length(obj.Nu0)).*obj.Nu0(1)];
            end
            nuMat = repmat(obj.Nu0, 2, 1);
            nuMat(2,:) = nuMat(2,:) .* 0.1;
%             nuMat = nuMat .* 0.1;
            obj.Nu0 = reshape(nuMat, [1, obj.nClasses*2]);
            
            obj.setAlphaPrior(bccSettings.AlphaDiff);
        end
                
        function initAlpha(obj, Alpha)
        
            if isempty(obj.Alpha0)
                obj.Alpha0 = zeros(obj.nClasses.*2, obj.nScores, obj.nAgents);

                obj.Alpha0(:, :, :) = 1 ./ obj.nScores;      
                
                %sequential version needs stronger priors to avoid assuming too much from the early data points
                if obj.sequential
                    obj.Alpha0 = Alpha + 2;
                end
            else
                if size(obj.Alpha0, 3) < obj.nAgents
                    if obj.trustFinalAgent > 0
                        obj.trustedAlpha = reshape(repmat(obj.trustedAlpha', 2, 1), [obj.nScores, obj.nClasses*2, size(obj.trustedAlpha,3)])';
                    end
                    obj.setAlphaPrior(obj.Alpha0, true, obj.nAgents);
                end
            end       
        end
        
        function [ET, Tind] = initET(obj, C, nAssets, lnP)
            ET = sparse(obj.nClasses*2, nAssets);
            Tind = sparse(obj.nClasses*2, nAssets);
            
            if iscell(C)
                ET(:, C{2}) = 1./(obj.nClasses*2);            
            else
                ET = ET + 1/obj.nClasses;
            end
            for j=1:obj.nClasses
                ET(:, obj.targets==j) = 0;
                Tind(:, obj.targets==j) = 0;
                nu = obj.Nu0;
                ET(j*2-1, obj.targets==j) = exp(lnP(j*2-1));
                Tind(j*2-1, obj.targets==j) = exp(lnP(j*2-1));
                ET(j*2, obj.targets==j) = exp(lnP(j*2));
                Tind(j*2, obj.targets==j) = exp(lnP(j*2));                
            end    
        end      
        
        function lnP = initLnP(obj)
            nu = obj.Nu0;
            
            priorEP = nu ./ sum(nu);
            
            nuSum = nu;
            
            for j=1:obj.nClasses
                nu(j*2-1) = nu(j*2-1) + priorEP(j*2-1).*sum(obj.targets==j);
                nu(j*2) = nu(j*2) + priorEP(j*2).*sum(obj.targets==j);
                nuSum(j*2-1) = nu(j*2-1) + nu(j*2);
                nuSum(j*2) = nu(j*2-1) + nu(j*2);
            end            
            lnP = psi(nu) - psi(nuSum);         
        end        
        
        function Xi = updateApproxXi(obj, Count, ET)
            %rowTotals = repmat(sum(Count,3), [1, 1, obj.nAgents]);
            rowTotals = repmat(sum(Count,2), 1, obj.nScores);
%             obj.B = rowTotals./size(Count,1) ./ 100;
            %B = obj.B
            tau = (1+obj.B) ./ (rowTotals + obj.B);

            Xi = zeros(obj.nClasses*2, obj.nScores, size(Count,3));
            for j=1:obj.nClasses
              tau_j = tau(2*j-1:2*j,:,:);
              Count_j = Count(2*j-1:2*j,:,:);
              rowTotals_j = rowTotals(2*j-1:2*j,:,:);

              Xi(2*j,:,:) = (sum(tau_j.*Count_j,1) + obj.beta.*obj.gamma) ./ (sum(tau_j.*rowTotals_j,1) + obj.beta);
              Xi(2*j-1,:,:) = Xi(2*j,:,:);
            end
        end
        
        function Count = voteCounts(obj, C, T, agentIdx)
            testVoteCounts = zeros(obj.nClasses*2, obj.nScores, obj.nAgents);    
            Count = zeros(obj.nClasses*2, obj.nScores, obj.nAgents);    
                        
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
                                
                for j=1:obj.nClasses*2
                    for s=1:obj.nScores
                        Count(j, s, :) = (obj.Cmat==s) * T(j,:)';%sum(repmat(T(j,:), obj.nAgents, 1) .* obj.Cmat==s, 2);
                    end
                end
                obj.trainVoteCounts = Count;
            end
            
            for j=1:obj.nClasses*2
                try
                    for l=1:obj.nScores
                        if l==obj.settings.zeroScoreMap
                            continue;
                        end
                        Count(j, l, :) = (obj.Cmat==l) * T(j,:)';
                    end
                    if obj.settings.zeroScoreMap > 0
                        testVoteCounts(j, obj.settings.zeroScoreMap, :) = sum(T(j,obj.testIdxs)) - sum(testVoteCounts(j,:,:),2);                    
                        Count(j, obj.settings.zeroScoreMap, :) = testVoteCounts(j, obj.settings.zeroScoreMap, :) + obj.trainVoteCounts(j, obj.settings.zeroScoreMap, :);1
                    end

                catch me
                    display(['voteCounts fail in IBCC ' me.message]);
                end
            end
            
            Count = obj.mergeLowCounts(Count);
        end       
                
        %replace ET with double the rows for each difficulty level
        function Alpha = updateAlpha(obj, C, ET, agentIdx, OldPostAlpha, OldPostAlpha_t)
            if nargin < 5
                Count = obj.voteCounts(C, ET);
                Alpha = obj.Alpha0 + Count;
                %Alpha = Alpha + obj.B .* obj.updateApproxXi(Alpha,ET);
            else
                Count_a = obj.voteCounts(C, ET, agentIdx);
                Alpha = OldPostAlpha; %modified from tested version which appeared to be wrong
                Alpha(:,:,agentIdx) = obj.Alpha0(:,:,agentIdx) + Count_a;
            end                 
        end
        
        function [ELnPi, Alpha, EEta] = expectedLnPi(obj, C, ET, agentIdx, OldPostAlpha, OldPostAlpha_t)
            %specifying agentIdx means we only update for that agent from
            %the last score in C and use the given counts for others

            if nargin < 5
                Alpha = obj.updateAlpha(C, ET);
            else
                Alpha = obj.updateAlpha(C, ET, agentIdx, OldPostAlpha, OldPostAlpha_t);                  
            end
            
            %this is required because some indexes will have zeros in
            %otherwise which won't work with psi function. However, not
            %sure we need to set default values here?
            if isempty(obj.defaultLnPi)
                defLnPi = psi(obj.Alpha0(:,:,1));
                defNormTerm = psi(sum(obj.Alpha0(:,:,1), 2));
                
                obj.defaultLnPi = zeros(obj.nClasses*2, obj.nScores, obj.nAgents);
                obj.defaultNormTerm = zeros(obj.nClasses*2, obj.nScores, obj.nAgents);
                for j=1:obj.nClasses*2
                    for l=1:obj.nScores
                        obj.defaultLnPi(j,l,:) = defLnPi(j, l);
                        obj.defaultNormTerm(j,l,:) = defNormTerm(j);
                    end
                end
            end            
            if iscell(C)
                nonDefAgents = unique(C{1}); %non-default agents
            else
                nonDefAgents = 1:obj.nAgents;
            end

            ELnPi = obj.defaultLnPi - obj.defaultNormTerm;            
            try
                ELnPi(:, :, nonDefAgents) = psi(Alpha(:, :,nonDefAgents));
            catch e
                display(['psi error: ' e.messsage]);
            end
                        
            normTerm = psi(sum(Alpha(:,:,nonDefAgents), 2));
            normTerm = repmat(normTerm, 1, obj.nScores);
            ELnPi(:,:,nonDefAgents) = ELnPi(:,:,nonDefAgents) - normTerm;
            EEta = [];
        end        

        function [ElnP, PostNu] = expectedLnP(obj, T)
            
            ElnP = zeros(obj.nClasses, 1);
            
            totalN = sum(sum(T,1)~=0);
            
            normTerm = psi(sum(obj.Nu0) + totalN);
            Nu = obj.Nu0;
            for j=1:obj.nClasses*2
                try
                    Nu(j) = obj.Nu0(j) + sum(T(j, :));
                    ElnP(j) = psi(Nu(j)) - normTerm;
                catch me
                    display(['bad Nu or T values: ' me]);
                end
            end  
            
%             ElnPP = gammaln(sum(Nu', 1))-sum(gammaln(Nu'),1) + sum((Nu'-1).*ElnP, 1);
            PostNu = Nu;
        end        
        
        function [ET, pT, lnpCT] = expectedT(obj,  C, lnKappa, lnPi, nAssets)     
            
            if nargin < 5
                nAssets = length(obj.targets);
            end
            
            try
                %why do we use the zero score map like this? If we are
                %instantiating all the valid decisions that are given a
                %score of zero, why not just put them in the scoreset C to
                %start with? Adding obj.settings.zeroScoreMap to all scores at
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
            lnpCT = sparse([], [], [], nAssets, obj.nClasses*2, length(C{1}));
            pT = sparse([], [], [], nAssets, obj.nClasses*2, length(C{1})*obj.nClasses);

            for j=1:obj.nClasses*2
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

            expA = sparse(size(lnpCT,1), obj.nClasses*2);
            for j=1:obj.nClasses*2
                rescale = -min(lnpCT(:, :),[],2);

                expAj = exp(lnpCT(:,j)+ realIdxMat.*rescale);
                %stop any extreme values from causing Na
                expAj(expAj==Inf) = 1;
                expA(:,j) = expAj;
            end
            tooBig = find(sum(expA>10^100,2));
            expA(tooBig, :) = expA(tooBig, :) ./ repmat(max(expA(tooBig,:),[],2), 1, size(expA,2));
            expB = sum(expA,2);

            %using loops so we avoid having to repmat with large matrices
            for j=1:obj.nClasses*2
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
            
            ET = sparse(obj.nClasses*2, size(pT,2));
            ET(:, obj.testIdxs) = pT(obj.testIdxs, :)';
            
            for j=1:obj.nClasses
                pT(obj.trainIdxs,j*2-1) = pT(obj.trainIdxs,j*2-1) .* (obj.Tmat(j,obj.trainIdxs))';
                pT(obj.trainIdxs,j*2) = pT(obj.trainIdxs,j*2) .* (obj.Tmat(j,obj.trainIdxs))';
            end
            pT(obj.trainIdxs,:) = pT(obj.trainIdxs,:) ./ repmat(sum(pT(obj.trainIdxs,:),2), 1, obj.nClasses*2);
            ET(:, obj.trainIdxs) = pT(obj.trainIdxs, :)';
            if sum(sum(isnan(ET)))>0
                display('nan!');
            end
%             ET = pT;
%             ET(:, obj.targets~=0) = 0;
%             ET(sub2ind(size(ET), double(obj.targets(obj.targets~=0)), find(obj.targets~=0) )) = 1;
%             ET(1, obj.targets~=0) = 1-ET(2, obj.targets~=0);
        end        

        function combinedPost = combineCluster(obj, post, Tvec)
            
            obj.setTargets(Tvec);
            
            obj.nKnown = round(obj.nKnown);
            C = obj.prepareC(post);
            
            if iscell(C)
                obj.Cmat = sparse(C{1}, C{2}, C{3}, obj.nAgents, length(Tvec));
            else
                obj.Cmat = C;
            end
            %merge any base classifiers with too little information
%             C = obj.mergeRareBase(C);
                
            obj.initAlpha();             
            [ET, nIt, EPi, EAlpha] = obj.combineScoreSet(C);  
            obj.combinedPost = zeros(obj.nClasses, size(ET,2));
            for j=1:obj.nClasses
                obj.combinedPost(j,:) = ET(j*2-1,:) + ET(j*2,:);
            end
            
            if obj.nClasses==2
                combinedPost(1, :) = obj.combinedPost(2, :);
            else
                [maxVals, maxIdx] = max(ET, [], 1);
                combinedPost = obj.combinedPost;%ET(sub2ind(size(ET), maxIdx, 1:length(maxIdx))) + maxIdx - 1;
            end

            nItFile = sprintf('%s%s', obj.confSave, 'nIts_vb');
            if exist(nItFile, 'file');
                prevNIts = dlmread(nItFile);
                prevNIts = [prevNIts nIt];
            else
                prevNIts = nIt;
            end
%             dlmwrite(nItFile, prevNIts);
            
            display([num2str(nIt) ' iterations for IBCCVB. sequential=' num2str(obj.sequential)]);
            
%             obj.Alpha = obj.Alpha + EAlpha;
        end        
        
    end
end