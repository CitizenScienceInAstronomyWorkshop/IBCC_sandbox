classdef  CountDecisions < combiners.AbstractCombiner
    %SELECT_BEST Trivial "combiner" that simply determines the base learner
    %with lowest error over a window of earlier data points and uses its outputs.
    % Assumes there is online but maybe occasional access to target labels,
    % and uses the available ones from the most recent n test samples. If
    % there are less than min number of targets over the period n, the
    % method will look further back to find more labels.
    
    properties
        label = 'Decision Count';
        
    end
    properties (Constant)
        shortLabel = 'Decision Count';
    end
    
    methods
        function obj = CountDecisions(nAgents, K)
            obj@combiners.AbstractCombiner(nAgents, K);
            obj.normalised = 1;
        end
        
        function combinedPost = combineDecisions(obj, baseOutputs, clusters)
                        
            nSamples = size(baseOutputs, 2);
            combinedPost = zeros(obj.K, nSamples);
            
            if obj.K == 1
                totalWeights = sum( baseOutputs~=obj.noScore, 1);
                combinedPost(totalWeights>0) = totalWeights(totalWeights>0);
            else            
                for k=1:obj.K
                    membership = clusters(:,:)==k;
                    combinedPost(k, :) = sum(membership, 1);
                end
            end
            obj.combinedPost = combinedPost;
        end
    end
    
end

