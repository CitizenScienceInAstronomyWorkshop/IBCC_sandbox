classdef  IbccVbEmInit < combiners.bcc.IbccVb
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
        initer
        flipPi = false;
    end
 
    methods (Static)
        function sl = shortLabel
            sl = 'Ibcc-VB-EMinit';
        end   
    end
    
    methods       
        function obj = IbccVbEmInit(nAgents, K, targets, agents, nClasses, nScores, nu)
            obj@combiners.bcc.IbccVb(nAgents, K, targets, agents, nClasses, nScores, nu);
            obj.label = 'VB I.B.C.C. with EM init'; 
        end        
        function [lnPi, lnK, ET, Tind] = initVariables(obj,  C, nAssets)    
            
            if obj.debug
                display('IBCC-VB init');
            end
            
            if nargin < 4
                nAssets = length(obj.targets);
            end

            if sum(obj.Nu0) == 0
                obj.Nu0 = 0.5;
            end     

            obj.initer = combiners.bcc.IbccEm(obj.settings, obj.nAgents, obj.K, ...
                                obj.targets, obj.agents, obj.nClasses, obj.nScores);
            obj.initer.scoreSet = obj.scoreSet;
            if nargin > 4
                obj.initer.combineCluster(C, obj.targets, members);
            else
                obj.initer.combineCluster(C, obj.targets);
            end 
            
            lnK = obj.initer.lnK';
            [ET, Tind] = obj.initET(C, nAssets, lnK);            
            lnPi = obj.initer.lnPi;
        end  
    end
        
    
end

