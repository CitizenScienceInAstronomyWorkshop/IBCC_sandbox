classdef IbccVbInitMeds < combiners.bcc.IbccVb
    %IBCCVBINITMEDS Overwrites the init method.
    
    methods(Static)
        function sl = shortLabel
            sl = 'IbccVbInitMeds';
        end        
    end
    
    methods                
        function [Alpha, lnPi, lnP, ET, Tind, Count] = initVariables(obj, Tvector, C, Alpha, nAssets)
            Alpha = obj.initAlpha(Alpha);       
            
            if nargin < 5
                nAssets = length(Tvector);
            end
            
            ET = sparse(obj.nClasses, nAssets);
            Tind = sparse(obj.nClasses, nAssets);
            for n=1:length(Tvector)
                %leave the unknown ones as equal probability of each
                if Tvector(n) == 0
                    %is there any score for this data point? If not, ignore
                    %it and set to zero
                    if sum(C{2}==n)==0
                        ET(:, n) = 0;
                    else
                        ET(:, n) = 1 ./ obj.nClasses;
                    end
                else
                    %with continuous multi-class this goes wrong
                    ET(Tvector(n), n) = 1;
                    Tind(Tvector(n), n) = 1;
                end
            end

            lnP = log(obj.nu ./ sum(obj.nu));
            lnPi = log(ones(obj.nClasses, obj.nScores, obj.nAgents) ./ obj.nScores);
            Count = obj.voteCounts(C, ET);

        end      
  
    end
    
end

