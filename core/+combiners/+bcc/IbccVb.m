classdef  IbccVb < combiners.bcc.Ibcc
 %T should be in matrix format with rows corresponding to classes, and
 %entries giving the probability of a class
 
 %Full posterior distributions? At the moment we get only the expected T values at the end.
 %While these do encode uncertainty in the correct label, the full distr. would
 %tell us whether this uncertainty comes from uncertainty about the correct
 %model (i.e. we haven't observed points like this before), or whether
 %we are certain that it is hard to classify this data point given these
 %inputs (i.e. we have observed this pattern lots before and can't find a way to
 %differentiate points of different classes).
    
    properties
        %discountAlpha = 0.99;  
        defaultLnPi = [];
        defaultNormTerm = [];
        
        fold = 1;
                        
        sequential = false;       
        nBootstrap = 1; %how many samples to use as a batch to bootstrap the parameters before running sequentially.    
        
        psiMap = [];
    end
 
    methods (Static)
        function sl = shortLabel
            sl = 'IbccVb';
        end   
                
        function sl = shortLabelSeq
            sl = 'IbccVbSeq';
        end
    end
    
    methods       
        function obj = IbccVb(bccSettings, nAgents, K, targets, agents, nClasses, nScores)
            % nAgents - number of base classifiers/agents to combine
            % K - number of clusters of base classifiers (set to 1 if you
            % are not doing any clustering
            % targets - a vector containing the target class labels.
            % Training labels should have values from 1 to nClasses; test
            % labels should have value 0.
            % agents - not used in IBCC, can set to []
            % nClasses - number of target classes
            % nScores - number of output values from the base classifiers.
            % Base classifier outputs usually come from the same set of
            % values as the target labels, in which case nClasses and
            % nScores are the same
            % nu - hyperparameter Nu0; e.g. [10 10] in a two-class problem
            % Alpha - hyperparemeter Alpha0; a 3-D matrix; 1st dimension
            % corresponds to target class; 2nd dimension to base classifier
            % output score; 3rd dimension to ID of base classifier
            obj@combiners.bcc.Ibcc(bccSettings, nAgents, K, targets, agents, nClasses, nScores);
            obj.label = 'VB I.B.C.C.'; 
            obj.nonDet = false;
            obj.psiMap = zeros(1000000,1);
        end
        
        function [ET, Tind] = initET(obj, C, nAssets, ~)
            
            Tind = sparse(obj.nClasses, nAssets);            
            
            if ~isempty(obj.combinedPost)
                ET = obj.combinedPost;
            else
            
                ET = sparse(obj.nClasses, nAssets);
                if iscell(C)
                    ET(:, C{2}) = 1./obj.nClasses;            
                else
                    ET = ET + 1/obj.nClasses;
                end
                for j=1:obj.nClasses
                    ET(:, obj.targets==j) = 0;
                    ET(j, obj.targets==j) = 1;
                end 
                obj.combinedPost = ET;
            end
            
            for j=1:obj.nClasses
                Tind(:, obj.targets==j) = 0;
                Tind(j, obj.targets==j) = 1;
            end   
        end
        
        function lnP = initLnP(obj)
            nu = obj.Nu0;
            for j=1:obj.nClasses
                nu(j) = nu(j) + sum(obj.targets==j);
            end            
            lnP = psi(nu) - psi(sum(nu));         
        end
        
        function [lnPi, lnP, ET, Tind] = initVariables(obj,  C, nAssets)    
            
            if obj.debug
                display('IBCC-VB init');
            end
            
            obj.initAlpha();
            
            if nargin < 4
                nAssets = length(obj.targets);
            end

            if sum(obj.Nu0) == 0
                obj.Nu0 = 0.5;
            end     

            lnP = obj.initLnP();
            
            [ET, Tind] = obj.initET(C, nAssets, lnP);

            Cknown = C;
            if iscell(C)
                knownIdxs = ismember(C{2}, find(obj.targets>0));
                for i=1:length(Cknown)
                    Cknown{i} = Cknown{i}(knownIdxs);
                end
            end
            lnPi = obj.expectedLnPi(Cknown, Tind);
        end      
        
        function initAlpha(obj, Alpha)
            if isempty(obj.Alpha0)
                obj.Alpha0 = zeros(obj.nClasses, obj.nScores, obj.nAgents);

                obj.Alpha0(:, :, :) = 1 ./ obj.nScores;      
                
                %sequential version needs stronger priors to avoid assuming too much from the early data points
                if obj.sequential
                    obj.Alpha0 = Alpha + 2;
                end
            else
                if size(obj.Alpha0, 3) < obj.nAgents
                    obj.setAlphaPrior(obj.Alpha0, true, obj.nAgents);
                end
            end 
            
            %Initialise the posterior
            obj.Alpha = zeros(obj.nClasses, obj.nScores, obj.nAgents);                   
        end
        
        function updateProgressMon(obj, ET, nIt)
            if ~isempty(obj.progressMonitor)
                obj.progressMonitor.combinerIteration(...
                    ET, obj.getId(), nIt, obj.nKnown); 
            end
        end
                 
        function [lnPi, lnP, ET, nIt] = iterate(obj, ...
            nItPrev, lnPi, lnP, C, nAssets, ET, agentIdx)
            converged = false;
            cIt = 0; nIt = nItPrev;
            %L = -inf;
            diff = 0;
            
            while ~converged   
                if obj.debug
                    display(['IBCC-VB iterations: ' num2str(nIt) ', diff: ' num2str(diff) ]);
                end
                
                %', lower bound: ' num2str(L)]);
                % oldL = L;
                oldET = ET;
                nIt = nIt + 1;
                                
                obj.updateProgressMon(ET,nIt);
                
                [ET, ~] = obj.expectedT(C, lnP, lnPi, nAssets);
                
                if nargin < 8 || agentIdx==0
                    lnPi = obj.expectedLnPi(C, ET);
                else            
                    lnPi = obj.expectedLnPi(C, ET, agentIdx, obj.Alpha);
                end
                lnP = obj.expectedLnP(ET);
                agentIdx = 0; %only use the single-agent update once, since any changes to ET that result from the new data point will require a change to the other agents
                
                % L = obj.lowerBound(ET,C,lnPi,lnP,PostAlpha);             
                diff = sum(sum(abs(ET - oldET)));

                if abs(diff) < obj.settings.convThreshold
                    cIt = cIt + 1;
                    %display(['converged at ' num2str(L)]);
                else
                    cIt = 0;
                end
                if ((cIt > obj.settings.convIt  || nIt-nItPrev > obj.settings.maxIt) && obj.settings.fixedNIt==0) ...
                        || (nIt >= obj.settings.fixedNIt && obj.settings.fixedNIt > 0)
                    converged = true;
                    % display(['IBCC-VB: final diff: ' num2str(diff)]);% ', ' num2str(L)]);
                end                
            end    
        end
                      
        function [ET, nIt] = combineScoreSet(obj, C)
                        
            if iscell(C)
                nSamples = length(C{2});
            else
                nSamples = 1;
            end
            nAssets = length(obj.targets);
            nIt = 0;
            %display(['IBCCVB: number of samples to process: ' num2str(nSamples)]);
            
            if obj.sequential
                firstIt = obj.nBootstrap;
                if firstIt > nSamples
                    firstIt = nSamples;
                end
            else
                firstIt = nSamples;                            
            end

            for i=firstIt:nSamples
            %display(['sample no: ' num2str(i) ', no. VB iterations: ' num2str(nIt)]);
                
                C_i = C;
                if iscell(C)
                    C_i{1} = C_i{1}(1:i);
                    C_i{2} = C_i{2}(1:i);
                    C_i{3} = C_i{3}(1:i);
                 
                    if i==0
                        currentAgent = 0;
                    else
                        currentAgent = C_i{1}(i); 
                        %if currentAgent is new and the sample is new we can
                        %simplify updates - this will not apply when we have to
                        %process multiple data points at once
                        if sum(C_i{1}==currentAgent)>1 || sum(C_i{2}==C_i{2}(i))>1 || i==firstIt
                            currentAgent = 0;
                        end
                    end
                end
                   
                if i==firstIt
                    [lnPi, lnP, ET] = obj.initVariables(C_i, nAssets);
                    currentAgent = 0;
                    obj.updateAlpha(C, ET);
                end
                %display('IBCC-VB begin iterator');                
                [lnPi, lnP, ET, nIt] = obj.iterate(nIt, lnPi, lnP, ...
                    C_i, nAssets, ET, currentAgent);                
            end

            EPi = exp(lnPi);
            obj.lnPi = lnPi;

            if obj.settings.debug
                obj.printPiEvolution(obj.Alpha, C);
            end
            
            obj.lnK = lnP;
            obj.Pi = EPi;
            obj.C = C;
        end
        
        function printPiEvolution(~, C)
            classifiers = unique(C{1});
            for k=classifiers'
                display(['Alpha for classifier ' num2str(k)]);
                obj.Alpha(:,:,k)
            end            
        end
        
        function updateAlpha(obj, C, ET, agentIdx, OldPostAlpha, ~)
            if nargin < 5
                Count = obj.voteCounts(C, ET);
                obj.Alpha = obj.Alpha0 + Count;              
            else
                Count_a = obj.voteCounts(C, ET, agentIdx);
                obj.Alpha = OldPostAlpha; %modified from tested version which appeared to be wrong
                obj.Alpha(:,:,agentIdx) = obj.Alpha0(:,:,agentIdx) + Count_a;
            end                 
        end
        
        function ELnPi = expectedLnPi(obj, C, ET, agentIdx, OldPostAlpha, OldPostAlpha_t)
            %specifying agentIdx means we only update for that agent from
            %the last score in C and use the given counts for others

            if nargin < 5
                obj.updateAlpha(C, ET);
            else
                obj.updateAlpha(C, ET, agentIdx, OldPostAlpha, OldPostAlpha_t);                  
            end
            
            %this is required because some indexes will have zeros in
            %otherwise which won't work with psi function. However, not
            %sure we need to set default values here?
            if isempty(obj.defaultLnPi)
                defLnPi = psi(obj.Alpha0(:,:,1));
                defNormTerm = psi(sum(obj.Alpha0(:,:,1), 2));
                
                obj.defaultLnPi = zeros(obj.nClasses, obj.nScores, obj.nAgents);
                obj.defaultNormTerm = zeros(obj.nClasses, obj.nScores, obj.nAgents);
                for j=1:obj.nClasses
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
            Alf = obj.Alpha(:,:,nonDefAgents);
            
            try
                mapIdxs = ceil(10000.*Alf);
                
                empties = mapIdxs>numel(obj.psiMap);
                obj.psiMap(mapIdxs(empties)) = psi(Alf(empties));

                empties = obj.psiMap(mapIdxs)==0;
                obj.psiMap(mapIdxs(empties)) = psi(Alf(empties));
                
                ELnPi(:, :, nonDefAgents) = obj.psiMap(mapIdxs);
                                
            catch e
                display(['psi error: ' e.message]);
            end
                      
            Alf = sum(Alf,2);
            mapIdxs = ceil(10000.*Alf);
            empties = mapIdxs>numel(obj.psiMap);
            obj.psiMap(mapIdxs(empties)) = psi(Alf(empties));
            empties = obj.psiMap(mapIdxs)==0;
            obj.psiMap(mapIdxs(empties)) = psi(Alf(empties));
            
            normTerm = obj.psiMap(mapIdxs);
            normTerm = repmat(normTerm, 1, obj.nScores);
            ELnPi(:,:,nonDefAgents) = ELnPi(:,:,nonDefAgents) - normTerm;
        end
        
        function ElnP = expectedLnP(obj, T)
            
            ElnP = zeros(obj.nClasses, 1);
            
            totalN = sum(sum(T,1)~=0);
            
            normTerm = psi(sum(obj.Nu0) + totalN);
            Nu = obj.Nu0;
            for j=1:obj.nClasses
                try
                    Nu(j) = obj.Nu0(j) + sum(T(j, :));
                    ElnP(j) = psi(Nu(j)) - normTerm;
                catch me
                    display(['bad Nu or T values: ' me]);
                end
            end  
            
            obj.Nu = Nu;
        end
        
        function [combinedPost, agentRatings] = combineDecisions(obj, C, ~) 
            % C - a cell array with 3 cells; the first contains
            % an ordered vector of base classifier IDs; the second contains a vector
            % of object IDs for the objects to be classified; the third
            % contains the base classifier output scores. 
            % clusters - set to [] if not performing clustering            
            obj.setTargets(obj.targets);
            obj.nKnown = round(obj.nKnown);
            C = obj.prepareC(C);
                                        
            [ET, nIt] = obj.combineScoreSet(C);
            
            obj.combinedPost = ET;
            if obj.nClasses==2
                combinedPost(1, :) = ET(2, :);
            else
                combinedPost = ET;
            end

            display([num2str(nIt) ' iterations for IBCCVB. sequential=' num2str(obj.sequential)]);
            
            agentRatings = obj.Alpha;            
        end
                      
        function [ET] = flipCounts(obj, ET, pT, ~, C, ~, ~)
            invET = ET;
            invET(1,obj.testIdxs) = ET(2,obj.testIdxs);
            invET(2,obj.testIdxs) = ET(1,obj.testIdxs);          
            invLnPi = obj.expectedLnPi(C, invET);
            invLnP = obj.expectedLnP(invET);

            [invET, invPT] = obj.expectedT(C, invLnP, invLnPi, size(ET,2));
            
            %calculate probabilities of ET and invET for weighting; use
            %joint distribution with probability of lnPi and lnP
            %not quite right as we should be using expectation of pi
            %and p rather than of log pi and log p.
                        
            %need to rescale
            lnpInvT = sum(log(sum(invPT.*invET,1)));
            lnpT = sum(log(sum(pT.*ET,1))) - lnpInvT;
            if isinf(exp(lnpT))
                pPeak1 = 1;
                pPeak2 = 0;
            elseif isinf(exp(lnpInvT))
                pPeak1 = 0;
                pPeak2 = 1;
            else
                pPeak1 = exp(lnpT) / (exp(lnpT)+1);
                pPeak2 = 1 / (exp(lnpT)+1);
            end

            ET = ET.*0.5 + invET.*0.5;
        end   
        
        %lower bound decreases sometimes occur: this only happens when both
        %p (or kappa) and pi are being updated - when one is fixed we don't
        %have a problem.
        function [L, EEnergy, H] = lowerBound(obj, T, C, lnPi, lnP, PostAlpha)

            J = size(lnPi,1);
            
            nResponses = length(C{1});
            ElnPC = zeros(J, nResponses);            
            for j=1:J
                scores = zeros(length(C{3}), 1) + obj.settings.zeroScoreMap;
                scores(C{3}~=0) = C{3}(C{3}~=0);
                idxs = sub2ind([size(lnPi,2) size(lnPi,3)], scores, C{1});                                
                ElnPC(j, :) = T(j,C{2}).*lnPi(j, idxs);
            end
            ElnPC = sum(sum(ElnPC));
            
            if size(lnP, 1) > size(lnP, 2)
                lnP = lnP';
            end
            TargetCounts = sum(T, 2);
            ElnPT = sum(TargetCounts .* lnP', 1);
            
            ElnPPi = gammaln(sum(obj.Alpha0, 2))-sum(gammaln(obj.Alpha0),2) + sum((obj.Alpha0-1).*lnPi, 2);
            ElnPPi = sum(sum(ElnPPi, 1), 3);
            
            PriorNu = obj.Nu0;
            ElnPP = gammaln(sum(PriorNu, 2))-sum(gammaln(PriorNu),2) + sum((PriorNu-1).*lnP, 2);
            
            EEnergy = ElnPC + ElnPT + ElnPPi + ElnPP;
        
            nAssets = size(T,2);
            nonZeroScores = C{3}~=obj.noScore;
            nValid = sum(nonZeroScores);
            validDecs = C{3}(nonZeroScores);
            validN = C{2}(nonZeroScores);
            validAgents = C{1}(nonZeroScores);
            indx = sub2ind([obj.nScores obj.nAgents], ...
                            validDecs, ...
                            validAgents);     
            lnpT = sparse(J, nAssets);
            for j=1:J
                lnpT(j,:) = sparse(ones(1,nValid), validN', lnPi(j, indx), 1, nAssets) + lnP(j);                
            end         
            lnpT(:, obj.targets~=0) = -inf;
            lnpT(sub2ind(size(T), obj.targets(obj.targets~=0), find(obj.targets~=0))) = 0;
            ElnQT = sum(sum(T(T~=0) .* lnpT(T~=0)));
            ElnQTTrue = sum(sum(T(T~=0) .* log(T(T~=0))));

            %ignore default values!
            ElnQPiTrue = gammaln(sum(PostAlpha, 2))-sum(gammaln(PostAlpha),2) + sum((PostAlpha-1).*lnPi, 2);
            ElnQPiTrue = sum(sum(ElnQPiTrue));
            ElnQPi = gammaln(sum(obj.Alpha0, 2))-sum(gammaln(obj.Alpha0),2) + sum((PostAlpha-1).*lnPi, 2);
            ElnQPi = sum(ElnQPi, 1);
            ElnQPi = sum(ElnQPi, 3);            
            
            TargetCounts = sum(T, 2);
            PostNu = obj.Nu0 + TargetCounts';
            if size(lnP, 1) > size(lnP, 2)
                lnP = lnP';
            end
            PriorNu = obj.Nu0;
            ElnQPTrue = gammaln(sum(PostNu, 2))-sum(gammaln(PostNu),2) + sum((PostNu-1).*lnP, 2);
            ElnQPTrue = sum(sum(ElnQPTrue));
            ElnQP = gammaln(sum(PriorNu, 2))-sum(gammaln(PriorNu),2) + sum((PostNu-1).*lnP, 2);
            
            H = - ElnQTTrue - ElnQPiTrue - ElnQPTrue;
            
            L = EEnergy + H;
        end     
    end
end
