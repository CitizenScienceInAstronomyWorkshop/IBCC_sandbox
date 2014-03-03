classdef  PrecWeightedAveraging < combiners.weighted.WeightedSum
    %PRECWEIGHTEDAVERAGING Averages the posteriors using weights derived
    %from precision of logodds. *This doesn't really make sense as there is
    %no reason to use precision of logodds to weight posteriors. Posteriors
    %are not Gaussians, but should be betas, so cannot find maximum 
    %likelihood by precision-weighted addition.*
    
    properties
        
    end
    properties (Constant)
        subShortLabel = 'precWeightedAverage';
    end
    
    methods
        function obj = PrecWeightedAveraging(nAgents, K, targets, agents, windowSize)
            obj@combiners.weighted.WeightedSum(nAgents, K, targets, agents, windowSize);
            obj.label = 'Model Averaging (precision weights)';
        end
        
        function createWeightCalculator(obj)
            obj.weightCalculator = combiners.weighted.weights.PrecisionWeights(obj.agents);
        end
    end
end

