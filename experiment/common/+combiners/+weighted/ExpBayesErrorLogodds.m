classdef ExpBayesErrorLogodds < combiners.weighted.NaiveBayesLogodds
    %EXPBAYESERRORLOGODDS Like NaiveBayesLogodds, but we weight the logodds by the
    %Expected Bayes Error rather than true precision.
        
    properties (Constant)
        subShortLabel = 'expectedBayesErrorWeightedLogodds';
    end
    
    methods
        function obj = ExpBayesErrorLogodds(nAgents, nClusters, labels, agents, windowSize)
            obj@combiners.weighted.NaiveBayesLogodds(nAgents, nClusters, labels, agents, windowSize);
        end
        
        function createWeightCalculator(obj)
            obj.weightCalculator = combiners.weighted.weights.ExpectedBayesErrorWeights(length(obj.agents));
        end
    end
    
end

