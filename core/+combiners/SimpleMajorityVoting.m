classdef SimpleMajorityVoting < combiners.AbstractCombiner
    %TRUE_BEST_AGENT This class selects that base learner that produced
    %least error for each data point. It is not a real combination or
    %selection method as it uses the true target label after each data
    %point, so is just used to compare with classifier combination
    %techniques.
     
    properties
        label = 'Maj. Voting'; 
        roundCombination = true;
        voteThreshold = 0.5; %boundary between for and against votes
        minScore = 0;
        maxScore = 1;
    end
    properties (Constant)
        shortLabel = 'majVote';
        shortLabelUnrounded = 'maj vote unrounded';
    end
    
    methods
        function obj = SimpleMajorityVoting(nAgents, K, nClasses)
            obj@combiners.AbstractCombiner(nAgents, K, nClasses);
            obj.normalised = true;
            obj.nClasses = nClasses;
        end
        
        function [combinedPost, agentRatings] = combineDecisions(obj, baseOutputs, clusters)
                         
            agentRatings = []; %we don't rate agents in this combiner           
            
            baseOutputs = round(baseOutputs);
            %number of voters - may vary between samples
            %nVoters = votes~=obj.noScore;
            
            %vote threshold is always 0 
            validIdxs = baseOutputs~=obj.noScore;
            if length(validIdxs) > 0.5*size(baseOutputs,1)*size(baseOutputs,2)
                baseOutputs = (baseOutputs-obj.minScore)/(obj.maxScore-obj.minScore) - 0.5;                
            else
                baseOutputs = (baseOutputs-validIdxs.*obj.minScore)./(obj.maxScore-obj.minScore) - 0.5.*validIdxs;
            end
%             
%             baseOutputs(baseOutputs==1) = 1;
%             baseOutputs(baseOutputs>1) = 2;
%             baseOutputs(baseOutputs<0) = -1;
            
            for k=1:obj.K  
                
                if obj.nClasses==2
                    if obj.K > 1
                        clusterMembers = (clusters==k);

                        voteSum = sum(baseOutputs .* clusterMembers, 1);
                    else
                        voteSum = sum(baseOutputs, 1);
                    end

                    %voteSum = voteSum ./ nVoters;
                    output = zeros(size(voteSum)) + 0.5;
                    output(voteSum>0) = 1;
                    output(voteSum<0) = 0;

                    obj.combinedPost(k, :) = output;
                else
                    N = hist(baseOutputs);
                    [maxVal obj.combinedPost(k, :)] = max(N, [], 1);
                    obj.combinedPost(k, :) = obj.combinedPost(k, :);
                end
            end
            combinedPost = obj.combinedPost;
            
            if obj.nClasses > 2
                combinedPost = zeros(obj.nClasses, size(combinedPost,2));
                for j=1:obj.nClasses
                    combinedPost(j,:) = obj.combinedPost==j;
                end
                
                obj.combinedPost = combinedPost;
            end
            
            if ~obj.roundCombination
                combinedPost = totalVotes1 ./ size(baseOutputs, 1);
            end        
        end        
    end
    
end
