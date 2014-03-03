classdef IbccVbAux < combiners.bcc.IbccVb
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        maxAlphaIt = 50;
    end
    
    methods (Static)
        function sl = shortLabel
            sl = 'IbccVbAux';
        end   
    end    
    
    methods
        function obj = IbccVbAux(nAgents, K, targets, agents, nClasses, nScores, ...
                lambdaMag, lambdaSym, nu)
            obj@combiners.bcc.IbccVb(nAgents, K, targets, agents, nClasses, nScores, ...
                lambdaMag, lambdaSym, nu);

            obj.label = 'VB-Aux I.B.C.C.'; 
        end
                
        function EAlpha = expectedAlpha(obj, Count, EAlpha)
            %in this version we do not use Pi to update Alpha as this
            %involves a complex conditional. For an approximation to the
            %full model of hyperparameter probabilities using auxiliary 
            %variables, see IbccVbAux.
            
            %seed Geometric mean of Alpha with its current value
            GAlpha = EAlpha;
            
            nIt = 0;
            
            OldEAlpha = EAlpha;
            
            while nIt < obj.maxAlphaIt % && OldEAlpha ~= EAlpha
                nIt = nIt + 1;

                EBeta = EAlpha;
                CountSum = Count;
                for k=1:obj.nAgents
                    for j=1:obj.nClasses
                        for i=1:obj.nScores
                            EBeta(j, i, k) = sum(EAlpha(j, :, k), 2) - EAlpha(j, i, k);
                            CountSum(j, i, k) = sum(Count(j, :, k), 2);
                        end
                    end
                end

                EM = GAlpha .* (psi(GAlpha+Count)-psi(GAlpha));
                EM(find(GAlpha==0)) = 0;
                EM(find(isinf(psi(GAlpha)))) = 0;
                %shouldn't this be log of expecation rather than expecation
                %of log?
                ELnOmega = log((EAlpha+EBeta) ./ (EAlpha+EBeta+CountSum));%psi(EAlpha+EBeta) - psi(EAlpha+EBeta+CountSum);

                GAlpha = exp(psi(EM)-log(obj.lambda-ELnOmega));

                OldEAlpha = EAlpha;
                EAlpha = (EM+1) ./ (obj.lambda-ELnOmega);
            end
        end    
        
%         function [ELnPi, Count] = expectedLnPi(obj, C, T, AlphaPrior)
%             PostAlpha = zeros(obj.nClasses, obj.nScores, obj.nAgents);
%             Count = zeros(obj.nClasses, obj.nScores, obj.nAgents);
%             ELnPi = zeros(obj.nClasses, obj.nScores, obj.nAgents);
%             
%             dupMatrix = ones(obj.nAgents, 1);
%             
%             for j=1:obj.nClasses
%                 
%                 Tj = dupMatrix * T(j, :);
%                 
%                 for l=1:obj.nScores
%                     Countjl = sum( (C==l) .* ( Tj ), 2);
%                     Alphajl = AlphaPrior(j, l, :);
%             
%                     Count(j, l, :) = Countjl;
%                     PostAlpha(j, l, :) = Count(j, l, :) + Alphajl;
%                 end
%                 
%                 normTerm = psi(sum(PostAlpha(j, :, :), 2));
% 
%                 for l=1:obj.nScores
%                     ELnPi(j, l, :) = psi(PostAlpha(j, l, :)) - normTerm;
%                 end
%             end
%         end        
    end
    
end

