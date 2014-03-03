classdef CIbccVb < combiners.bcc.IbccVb
    properties
        useNegFeatures = true; %count non-occurrences of features as a ngeative response
        useLikelihood = false;
        
        weight = 1;
        negC = [];
    end
    
    methods (Static)
        function sl = shortLabel
            sl = 'CIbccVb';
        end  
    end
    
    methods
        function obj = CIbccVb(bccSettings, nAgents, K, targets, agents, nClasses, nScores)
            obj@combiners.bcc.IbccVb(bccSettings, nAgents, K, targets, agents, nClasses, nScores);

            obj.label = 'VB C.I.B.C.C.'; 
            if nScores ~= 2
                display('continuous IBCC does not work with more than 2 scores at the moment');
            end
        end
        
        function mappedC = prepareC(obj, post)
            mappedC = post;            
        end        
        
        function Count = voteCounts(obj, C, T)
            Count = zeros(obj.nClasses, obj.nScores, obj.nAgents);
%             dupMatrix = ones(obj.nAgents, 1);
%             C = C.*obj.weight;
            
            Cpos = C; %Cpos(C<0) = 0;
%             CnoAnswer = C==-1;
%             Cneg = 1-C; %Cneg(C<0) = 0;

%             Cneg = C; Cneg(C<0) = 1; Cneg(C>weight) = weight; %prevent negative counts
%             Cneg = 1-Cneg;
            
            for j=1:obj.nClasses
                
                Tj = T(j, :);
                totalCounts = sum(Tj);
                Count(j, 2, :) = Cpos*Tj';
%                 if ~isempty(find(CnoAnswer, 1))
%                     Count(j, 1, :) = - CnoAnswer*Tj';
%                 end
                Count(j, 1, :) = Count(j, 1, :) + totalCounts - Count(j, 2, :);
            end
        end     

        function [ET, pT, logJoint] = expectedT(obj, C, lnP, lnPi, ~)
            logJoint = zeros(obj.nClasses, size(C, 2));
                
%             C = obj.weight.*C;
            Cpos = C; %Cpos(C<0) = 0;
            if obj.useNegFeatures && isempty(obj.negC)
                obj.negC = obj.weight-Cpos;
            end
            
            for j=1:obj.nClasses      
                %Turn -1s into zeros when totting up the positive examples
                %and -1s into 1s when totting up the negative examples

                if obj.useNegFeatures
                    lnPCT = reshape(lnPi(j,2,:), 1, obj.nAgents)*Cpos + ...
                        reshape(lnPi(j,1,:), 1, obj.nAgents)*(obj.negC);
                else
                    lnPCT = reshape(lnPi(j,2,:), 1, obj.nAgents)*Cpos;
                end
                
                if obj.useLikelihood
                    logJoint(j, :) = lnPCT;
                else

                    logJoint(j, :) = lnP(j) + lnPCT;
                end
            end
            rescale = repmat(-min(logJoint,[],1), obj.nClasses, 1);
            pT = exp(logJoint+rescale);
            normTerm = ones(obj.nClasses, 1)*sum(pT, 1);
            pT = pT ./ normTerm;
                        
            ET = obj.Tmat;
            ET(:, obj.testIdxs) = pT(:, obj.testIdxs);
        end        
    end
    
end

