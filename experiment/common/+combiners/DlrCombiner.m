classdef  DlrCombiner < combiners.AbstractCombiner
    %SELECT_BEST Trivial "combiner" that simply determines the base learner
    %with lowest error over a window of earlier data points and uses its outputs.
    % Assumes there is online but maybe occasional access to target labels,
    % and uses the available ones from the most recent n test samples. If
    % there are less than min number of targets over the period n, the
    % method will look further back to find more labels.
    
    properties
        label = 'DLR';
                roundCombination = false;
    end
    properties (Constant)
        shortLabel = 'dlr';
        shortLabelRounded = 'dlr rounded';

    end
    
    methods
        function obj = DlrCombiner(nAgents, K, targets, dlrPath)
            obj@combiners.AbstractCombiner(nAgents, K, targets);
            addpath(dlrPath);
        end
        
        function [combinedPost, agentRatings] = combineDecisions(obj, baseOutputs, clusters)
               
            nClasses = length(unique(obj.targets)) - 1;
               
            if nClasses > 2
                combinedPost = zeros(nClasses, length(baseOutputs));
            else
                nClasses = nClasses-1;
                combinedPost = zeros(1, length(baseOutputs));
            end
            
            dlrs = {nClasses};
            
            T = obj.targets;
            if length(T) < size(baseOutputs, 2)
                T = [T -1*ones(1, size(baseOutputs, 2)-length(T))]';
            else
                T = T(1:size(baseOutputs,2));
            end            
            
            %get the relevant inputs
            for j=1:nClasses
                Tj = sparse(T==j);
                dlrs{j} = dlr(baseOutputs', Tj);    
            end
            
            for j=1:nClasses
                os = dlrs{j};

                combinedPost(j, :) = os.y;

%                 if length(os.a) ~= length(combinedLogodds)
%                     display('mismatch');
%                 end
            end

            if obj.roundCombination
                combinedPost = round(combinedPost);
            end
            
            agentRatings = os.w;
        end
    end
    
end

