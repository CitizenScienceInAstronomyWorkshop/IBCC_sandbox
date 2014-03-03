classdef TrueBestAgent < combiners.AbstractCombiner
    %TRUE_BEST_AGENT This class selects that base learner that produced
    %least error for each data point. It is not a real combination or
    %selection method as it uses the true target label after each data
    %point, so is just used to compare with classifier combination
    %techniques.
    %Warning: this is ANTI-CAUSAL!
    
    properties
        label = 'True Best Agent'; 
    end
    properties (Constant)
        shortLabel = 'trueBestAgent';
    end
    
    methods
        function obj = TrueBestAgent(targets, nAgents, K)
            obj@combiners.AbstractCombiner(nAgents, K, targets);
        end
        
        function combinedPost = combineDecisions(obj, baseOutputs, clusters)
            
            T = ones(obj.nAgents, 1) * obj.targets;
            
            differences = abs(T-baseOutputs);
                        
            obj.combinedPost = zeros(obj.K, length(baseOutputs));
            for k=1:obj.K
                
                if obj.K == 1
                    [min_errors best_idx] = min(differences, [], 1);
                else                    
                    best_idx = zeros(1, length(baseOutputs));
                    for n=1:length(baseOutputs)
                        col = differences(:,n);
                        cluster_members = col(clusters(:,n)==k);
                        min_error = min(cluster_members);
                        best = find(col==min_error);
                        %choose the first one from best if there are more
                        %than one with the same amount of error.
                        best_idx(1, n) = best(1);
                    end
                end
            
                %change the best indices from column numbers to indices in a
                %1-D array.
                row_start_inc = 0:length(baseOutputs)-1;
                row_start_inc = row_start_inc .* obj.nAgents;
                best_array_idx = best_idx + row_start_inc;
            
                obj.combinedPost(k, :) = baseOutputs(best_array_idx);
            end
            combinedPost = obj.combinedPost;
        end        
    end
    
end

