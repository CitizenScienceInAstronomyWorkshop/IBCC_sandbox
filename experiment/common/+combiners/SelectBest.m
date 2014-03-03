classdef SelectBest < combiners.AbstractCombiner
    %SELECT_BEST Trivial "combiner" that simply determines the base learner
    %with lowest error over a window of earlier data points and uses its outputs.
    % Assumes there is online but maybe occasional access to target labels,
    % and uses the available ones from the most recent n test samples. If
    % there are less than min number of targets over the period n, the
    % method will look further back to find more labels.
    
    properties
        window_size
        min_window
        label = 'Select Min Bayes Error';
        ignoreLastAgent = false; %if the last agent is actually the target values
    end
    properties (Constant)
        shortLabel = 'minBayesError';
    end
    
    methods
        function obj = SelectBest(targets, nAgents, window_size, min_window, K, targetsAsAgents)
            obj@combiners.AbstractCombiner(nAgents, K, targets);
            obj.window_size = window_size;
            obj.min_window = min_window;
            obj.ignoreLastAgent = targetsAsAgents;
        end
        
        function combinedPost = combineDecisions(obj, baseOutputs, clusters)
                        
            T = ones(obj.nAgents, 1) * obj.targets;
                        
            differences = abs(T-baseOutputs(:, 1:length(obj.targets)));
            
            obj.combinedPost = zeros(obj.K, length(baseOutputs));
            
            if obj.K > 1
                obj.selectByCluster(baseOutputs, clusters, differences);
            else
                obj.selectFromAll(baseOutputs, differences);
            end
            combinedPost = obj.combinedPost;
        end
        
        function selectByCluster(obj, baseOutputs, clusters, differences)
            best_idx = ones(obj.K, length(baseOutputs));
            
            for k=1:obj.K
                best_idx(k, 1) = find(clusters(:, 1)==k, 1);
                obj.combinedPost(k, 1) = baseOutputs(best_idx(k, 1), 1);
            end
            for i=2:length(baseOutputs)
 
                finish = i-1;

                start = finish-obj.window_size;
                if start < 1
                    start = 1;
                end
                
                real_targets = find(obj.targets(start:finish) > -1);
                    
                %Work out how many target values we can use to assess which
                %is the best agent to select.
                num_real_targets = length(real_targets);
                if num_real_targets < obj.min_window
                    
                    real_targets = find(obj.targets(1:finish) > -1);
                    num_real_targets = length(real_targets);
                    
                                       
                    start_real_targets = num_real_targets-obj.min_window;
                    if start_real_targets < 1 
                        start_real_targets = 1;
                    end
                    if num_real_targets > 0
                        real_targets = real_targets(start_real_targets:num_real_targets);
                       
                    end
                end
                
                if num_real_targets > 0                        
                    sum_diffs_window = sum(differences(:, real_targets), 2);

                    for k=1:obj.K

                        cluster_members = sum_diffs_window(clusters(:, i)==k, :);
                        min_combined_error = min(cluster_members);
                        best = find(sum_diffs_window==min_combined_error);
                        best_idx(k, i) = best(1);

                        obj.combinedPost(k, i) = baseOutputs(best_idx(k, i), i);
                    end

                else
                    %pick a representative from each cluster at random as
                    %we have no real targets to assess them by at this
                    %stage.
                    for k=1:obj.K
                        best_idx(k, i) = find(clusters(:, 1)==k, 1);
                        obj.combinedPost(k, i) = baseOutputs(best_idx(k, i), i);
                    end                        
 
                end
            end
        end
        
        function selectFromAll(obj, baseOutputs, differences)
            best_idx = ones(obj.K, length(baseOutputs));
            
            %this is meant to deal with lack of target labels - it's not
            %right yet so come back to this later
            %differences = differences(:, find(obj.targets > -1));
             
            %first result has no previous target labels to use. Pick at
            %random
            best_idx(1, 1) = 1;
            obj.combinedPost(1,1) = baseOutputs(1, 1);

            for i=2:length(baseOutputs)
 
                finish = i-1;

                start = finish-obj.window_size;
                if start < 1
                    start = 1;
                end
                
                real_targets = find(obj.targets(start:finish) > -1);
                    
                %Work out how many target values we can use to assess which
                %is the best agent to select.
                num_real_targets = length(real_targets);
                if num_real_targets < obj.min_window
                    
                    real_targets = find(obj.targets(1:finish) > -1);
                    num_real_targets = length(real_targets);
                    
                                       
                    start_real_targets = num_real_targets-obj.min_window;
                    if start_real_targets < 1 
                        start_real_targets = 1;
                    end
                    if num_real_targets > 0
                        real_targets = real_targets(start_real_targets:num_real_targets);
                       
                    end
                end
                
                if num_real_targets > 0
                    sum_diffs_window = sum(differences(:, real_targets), 2);

                    [combined_error, current_best] = min(sum_diffs_window);
                    best_idx(1, i) = current_best;
                    obj.combinedPost(1, i) = baseOutputs(current_best, i);
                else
                    %pick a representative from each cluster at random as
                    %we have no real targets to assess them by at this
                    %stage.
                    best_idx(1, i) = 1;
                    obj.combinedPost(1,i) = baseOutputs(1, i);
                end
            end
            
            %change the best indices from column numbers to indices in a
            %1-D array.
            
            
            %for k=1:obj.K
            %    best_array_idx = sub2ind(size(baseOutputs), best_idx(k,:), 1:length(baseOutputs));                
                %obj.combinedPost(k, :) = baseOutputs(best_array_idx);
                %obj.combinedLogodds(k, :) = base_logodds(best_array_idx);
            %end
            
        end
    end
    
end

