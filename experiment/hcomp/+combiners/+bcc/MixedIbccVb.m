classdef MixedIbccVb < combiners.bcc.DynIbccVb
    %DYNIBCCVB Dynamic IBCC-VB
    %   For base classifiers that change. Uses a dynamic
    %   dirichlet-multinomial model to determine the confusion matrices
    %   after each data point.
    
    properties
        Cfeat = []; %features in matrix format
        featureIdxs = []; %static, continuous-response features
        dynIdxs = []; %dynamic, discrete-response base classifiers
        
        featCombiner = [];
        
        logJoint = [];
    end
    
    methods(Static)
        function sl = shortLabel
            sl = 'MixIbccVb';
        end        
    end
    
    methods
        
        function obj = MixedIbccVb(bccSettings, nFeat, featSettings, features, nAgents, K, targets, agents, nClasses, nScores)
            obj@combiners.bcc.DynIbccVb(bccSettings, nAgents, K, targets, agents, nClasses, nScores);
            obj.label = 'Dynamic/Static mixture VB I.B.C.C.'; 
            
            obj.featureIdxs = 1:nFeat;
            obj.Cfeat = features;
            obj.dynIdxs = nFeat:nAgents;   
            
            obj.featCombiner = combiners.bcc.CIbccVb(featSettings, nFeat, 1, targets, [], nClasses, 2);
            obj.featCombiner.useNegFeatures = false;
            obj.featCombiner.useLikelihood = true;
        end        
         
        function [ElnPiBundle, AlphaBundle] = expectedLnPi(obj, C, ET, ~, ~)
            [ElnPiBundle, AlphaBundle] = expectedLnPi@combiners.bcc.DynIbccVb(obj, C, ET);
            [ElnPiFeat, AlphaFeat] = obj.featCombiner.expectedLnPi(obj.Cfeat, ET);
            ElnPiBundle{3} = ElnPiFeat;
            AlphaBundle{3} = AlphaFeat;
        end

        function [ET, pT] = expectedT(obj, C, lnP, lnPiBundle, nAssets)     
            
            if nargin < 6
                nAssets = length(obj.targets);
            end
            
            nonZeroScores = find(C{3}~=0);
            nValid = length(nonZeroScores);
            
            validDecs = C{3}(nonZeroScores); % c
            validN = C{2}(nonZeroScores); % index of data points

            indx = sub2ind([obj.nScores length(nonZeroScores)], ...
                            validDecs, ...
                            nonZeroScores);     
            lnPT = zeros(obj.nClasses, nAssets);
            lnPi = lnPiBundle{1};%Eta - log(1+exp(Eta));
            for j=1:obj.nClasses                
                lnPT(j,:) = sparse(ones(1,nValid), validN', lnPi(j, indx), 1, nAssets);
                if ~obj.settings.useLikelihood
                    lnPT(j,:) = lnPT(j,:) + lnP(j);
                end
            end
            
            [~,~, lnPFT] = obj.featCombiner.expectedT(obj.Cfeat, lnP, lnPiBundle{3}, nAssets);
            lnPT = lnPT + lnPFT;
            
            rescale = repmat(-min(lnPT), obj.nClasses, 1);
            expA = exp(lnPT+rescale);%exp(lnPT(:,realIdxs)+rescale);
            
            obj.logJoint = lnPT;

            expB = repmat(sum(expA,1), obj.nClasses, 1);
            
            %stop any extreme values from causing Na
            expB(expA==Inf) = 1;
            expA(expA==Inf) = 1;
            
            pT = expA./expB;
            ET = pT;            
            if length(obj.targets)<1
                return
            end
            ET = obj.Tmat;
            ET(:, obj.testIdxs) = pT(:, obj.testIdxs);
        end
        
        function ET = combineCluster(obj, post, Tvec)   
            obj.setTargets(Tvec);       
            obj.nKnown = round(obj.nKnown);
            C = obj.prepareC(post);
            obj.initAlpha(length(C{1}));            
            [ET, ~] = obj.combineScoreSet(C); 
            obj.combinedPost = ET;
        end  
        
        function initAlpha(obj, nDupes, Alpha)
            if isempty(obj.Alpha0)
                obj.Alpha0 = zeros(obj.nClasses, obj.nScores, nDupes);

                obj.Alpha0(:, :, :) = 1 ./ obj.nScores;      
                
                %sequential version needs stronger priors to avoid assuming too much from the early data points
                if obj.sequential
                    obj.Alpha0 = Alpha + 2;
                end
            else
                if size(obj.Alpha0, 3) < nDupes
                    obj.setAlphaPrior(obj.Alpha0, true, nDupes);
                end
            end             
        end       
            
        function [L, EEnergy, H] = lowerBound(obj, T, C, lnPiBundle, lnP, AlphaBundle)

        end
    end 
end

