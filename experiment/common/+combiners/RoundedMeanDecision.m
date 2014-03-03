classdef RoundedMeanDecision < combiners.MeanDecision
    %RoundedMeanDecision Round the mean decision to one of three answers: ...
    % 0, 1 or I don't know.
    
    properties (Constant)
        subShortLabel = 'roundedMean';
    end
    
    methods
        function obj = RoundedMeanDecision(nAgents, K)
            obj@combiners.MeanDecision(nAgents, K);
        end        
        
        function combinedPost = combineDecisions(obj, baseOutputs, clusters)
            
            meanPost = obj.combineDecisions@combiners.MeanDecision(baseOutputs);
            combinedPost = meanPost;
            combinedPost(meanPost<0) = -1; %min(combinedPost);
            combinedPost(meanPost>0) = 1.5;
            combinedPost(meanPost>1.5) = 3; %max(combinedPost);
            combinedPost(meanPost==0) = 1;
            %combinedPost(combinedPost>0&combinedPost<2) = ...
            %    (max(combinedPost)-min(combinedPost))/2+min(combinedPost);
        end
    end
    
end

