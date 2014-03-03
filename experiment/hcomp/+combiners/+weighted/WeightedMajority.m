classdef WeightedMajority < combiners.weighted.WeightedSum
    %WEIGHTEDMAJORITY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        postScores = true %whether the base classifiers provide posterior probs or not - if not we don't touch
        nVotesPerClassifier = 1;
        voteThreshold = 0.5;
    end
    
    properties (Constant)
        subShortLabel = 'weightedMajority';
        subShortLabelUnrounded = 'weightedVotes';
    end
    
    methods
        
        function obj = WeightedMajority(nAgents, K, targets, agents, windowSize)
            obj@combiners.weighted.WeightedSum(nAgents, K, targets, agents);
            if nargin > 4
                obj.windowSize = windowSize;
            end
            obj.label = 'Weighted Maj. Voting';
            obj.roundCombination = true;      
            obj.normalised = true;
            display('warning: WM may not work with positive vote values');
        end
        
        function scores = getScores(obj, scores)
            if obj.postScores
                scores = round(scores);
                scores = scores - obj.minScore - ((obj.maxScore-obj.minScore)/2);
            end
        end
        
        function combinedPost = combineDataPoints(obj, scores, W)
            %normalise so that smallest weight will be one
            scores = obj.getScores(scores);
            weightedVotes = sum(W.*scores, 1);
            
            %Weighted majority produces discreet outputs.
            if obj.roundCombination
                combinedPost = weightedVotes;
                combinedPost(weightedVotes>0) = 1;
                combinedPost(weightedVotes<0) = 0;
                combinedPost(weightedVotes==0) = 0.5;
            else
                %shift back to 0.5 as the mid value
                combinedPost = weightedVotes + 0.5;                
            end              
        end        
        
        function combinedPost_n = combinationStep(obj, post, invW, n, nEnd)
            if nargin < 5
                nEnd = n;
            end
            
            if n > size(invW, 2) || nEnd > n
                %assume that no-scores will be empty in a sparse matrix or
                %already 0
                post = post(:, n:nEnd);
                invW = invW(:, end);
            else
                %assume no-scores will be empty in a sparse matrix or
                %already 0
                post = post(:, n);
                invW = invW(:, n);
            end
                
            %avoid dividing votes by using weights rather than inverse
            %weights. Normalise so that smallest value will be one, largest
            %value will be unbounded
            nW = max(invW) ./ invW; % largest thing we divide votes by is 1

            post = obj.getScores(post);
            weightedVotes = nW' * post;

            %Weighted majority produces discreet outputs.
            if obj.roundCombination
                combinedPost_n = weightedVotes;
                combinedPost_n(weightedVotes>0) = 1;
                combinedPost_n(weightedVotes<0) = 0;
                combinedPost_n(weightedVotes==0) = 0.5;
            else
                %shift back to 0.5 as the mid value
                combinedPost_n = weightedVotes + 0.5;                
            end         
            obj.combinedPost = [obj.combinedPost combinedPost_n];
        end
       
    end
    
end

