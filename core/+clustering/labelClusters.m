function [ cluster_labels ] = labelClusters( clusters, nAgents, K, means, clusterData )
%labelClusters Labels cluster with lowest magnitude mean as 0, cluster
%with highest mag as 1. I.e. cluster 0 is for the most uncertain cluster,
%cluster 1 for the cluster with strongest belief.
    
    %{
    iteration_members = zeros(K, nAgents);
            
    for i=1:nAgents

        for k=1:K

            if clusters(i, 1) == k
                iteration_members(k, i) = 1;
            end
        end
    end
    %}
   
    
    %find lowest magnitude of mean
    %{
    lowest_mean = inf;
    lowest_k = -1;
    
    for k=1:K
        if abs(means(k,1)) < lowest_mean
            lowest_mean = abs(means(k, 1));
            lowest_k = k;
        end
    end
    %}
    
    if strcmp(clusterData, 'post')
        means = means - 0.5;
    end
    
    means = abs(means);
    [sorted_means sorted_idx] = sort( sum(means, 2) );
    
    cluster_labels = sorted_idx;%zeros(K, 1);
    
    %{
    for k=1:K
        if k == lowest_k
            cluster_labels(lowest_k, 1) = 1;
            %cluster_members(1,  :) = iteration_members(k, :);
        else
            cluster_labels(k, 1) = 2; %only works for two clusters atm.            
            %cluster_members(2,  :) = iteration_members(k, :);
        end
    end
    %}
    %cluster_labels
    %cluster_labels(1:K, 1) = 1:K;
      
    %cluster_members = iteration_members; 
    
    %{
    %label the clusters consistently in each iteration by labelling the
    %clusters the same as clusters of most similar size in previous steps
    if cluster_members == zeros(K, nAgents)
        cluster_members = iteration_members;
        
        cluster_labels(1:K, 1) = 1:K;
    else
        for k=1:K
            
            min_diff = nAgents;         
            label = 0;
            for l=1:K
                %what is the difference in size between the clusters from
                %the previous data point?
                diff = sum(abs(iteration_members(k, 1:nAgents) - cluster_members(l, 1:nAgents)));
                if diff == 0
                    %same size - let's assume this is same cluster, so
                    %cluster k from this data point == cluster l from
                    %previous data point.
                    label = l;
                    break
                elseif diff < min_diff
                    %record the difference in size - the one with smallest
                    %difference (or the first one if they have same
                    %difference) is assumed to match cluster k from current
                    %data point
                    label = l;
                    min_diff = diff;
                end
            end
            
            cluster_labels(k, 1) = label;
            
            cluster_members(label, :) = iteration_members(k, :);
        end
    end
    %}
end

