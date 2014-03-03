classdef  MeanDecision < combiners.AbstractCombiner
    %SELECT_BEST Trivial "combiner" that simply determines the base learner
    %with lowest error over a window of earlier data points and uses its outputs.
    % Assumes there is online but maybe occasional access to target labels,
    % and uses the available ones from the most recent n test samples. If
    % there are less than min number of targets over the period n, the
    % method will look further back to find more labels.
    
    properties
        label = 'Mean';
        
        maxScore = 1;
        minScore = 0;
    end
    properties (Constant)
        shortLabel = 'mean';
    end
    
    methods
        function obj = MeanDecision(nAgents, K, nClasses, maxScore, minScore)
            obj@combiners.AbstractCombiner(nAgents, K, [], nClasses);
            
            if nargin > 3
                obj.maxScore = maxScore;
            end
            
            if nargin > 4
                obj.minScore = minScore;
            end
        end
        
        function [combinedPost, agentRatings] = combineDecisions(obj, baseOutputs, clusters)
                        
            nSamples = size(baseOutputs, 2);
            combinedPost = zeros(obj.K, nSamples);
            
            if obj.K == 1
                totalWeights = sum( baseOutputs~=obj.noScore, 1);
                combinedPost(totalWeights>0) = ...
                    sum(baseOutputs(:,totalWeights>0), 1) ./ totalWeights(totalWeights>0);
            else            
                for k=1:obj.K
                    membership = clusters(:,:)==k;
                    combinedPost(k, :) = sum( baseOutputs .* membership ) ./ sum(membership, 1);
                end
            end
            
            combinedPost = (combinedPost-obj.minScore) ./ (obj.maxScore-obj.minScore);
            
            if obj.nClasses > 2
                
                obj.combinedPost =  zeros(obj.nClasses, size(combinedPost,2));
                for j=1:obj.nClasses
                    obj.combinedPost(j,:) = 1 - abs(obj.combinedPost(j,:)-j);
                    roundIdxs = obj.combinedPost(j, :)<0;
                    obj.combinedPost(j,roundIdxs) = 0;
                end
                combinedPost = obj.combinedPost;
            else
                obj.combinedPost = combinedPost;
            end
            
            agentRatings = []; %we don't rate agents in this combiner
        end
    end
    
end

