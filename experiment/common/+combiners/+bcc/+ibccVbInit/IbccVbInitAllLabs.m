classdef IbccVbInitAllLabs < combiners.bcc.IbccVb
    %IBCCVBINITALLLABS Summary of this class goes here
    %   Detailed explanation goes here
    
    methods (Static)
        function sl = shortLabel
            sl = 'IbccVbInitAllLabs';
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

            if sum(obj.nu) == 0
                lnP = [log(0.5) log(0.5)];
            else
                lnP = psi(obj.nu) - psi(sum(obj.nu)); 
            end            

            [lnPi, Count] = obj.expectedLnPi(C, ET, Alpha, []);
        end           
        
    end
    
end

