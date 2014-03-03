classdef  NaiveBayesLogodds < combiners.weighted.WeightedCombiner
    
    properties
        label = 'N.B. LogOP'; %; base learners give indp. noisy estimates of a';
        X %input data to learners
        roundCombination = false;
    end   
    properties (Constant)
        shortLabel = 'naiveBayes';
        precalcShortLabel = 'naiveBayesPrec';
        precalcShortLabelRounded = 'nb rounded';
    end
    
    methods
        function obj = NaiveBayesLogodds(nAgents, K, targets, agents, windowSize, X)
            obj@combiners.weighted.WeightedCombiner(nAgents, K, targets, agents, windowSize); 
            if nargin > 5
                obj.X = X; %optional input data if we need to calculate precisions
            else
                obj.X = [];
            end
            
        end
        
        function createWeightCalculator(obj, ~, ~)
%             if isempty(obj.X)
%                 obj.weightCalculator = combiners.weighted.weights.PrecisionWeights(obj.agents, members);
%             else
                obj.weightCalculator = combiners.weighted.weights.LogOpWeights(obj.agents);
%             end
        end
        
        function combinedPost = combinationStep(obj, post, W)
            
            logodds = obj.weightCalculator.getLogodds(post);
            logodds = logodds - obj.weightCalculator.bias;
            
            %here we combine by logodds only. We find the naive Bayes
            %estimate for the logodds value.
            combinedLogodds = sum(repmat(W, 1, size(logodds,2)) .* logodds, 1);

            %apply the sigmoid function
            ebx = exp(-combinedLogodds);
            combinedPost = 1.0 ./ (ebx + 1.0);   
                
            obj.combinedPost = combinedPost;           
            
            if obj.roundCombination
                obj.combinedPost = round(obj.combinedPost);
                combinedPost = obj.combinedPost;
            end        
        end
    end
    
end

