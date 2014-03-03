classdef agentClustererKmeans < handle

    properties
        %cluster_members
        cluster_means
        cluster_errors
        cluster_sizes
        cluster_confidences
        
        agentClusterMembership
        
        sensor_cluster_credit
        sensor_cluster_weight
        sensor_cluster_acc_credit   
        
        weightClusterInput
        
        N
        A
        K
        L
        
        classifications
        
        clusterData
    end
    
    methods
        function obj = agentClustererKmeans(expSettings)
            
            import clustering.*
            
            obj.N = expSettings.nSamples;
            N = obj.N;
            
            obj.A = expSettings.nAgents;
            A = obj.A;
            
            obj.K = expSettings.K;
            K = obj.K;
            
            obj.L = expSettings.seqLength;
            
            obj.weightClusterInput = expSettings.weightClusterInput;
            
            %obj.cluster_members = zeros(K, A);
            obj.agentClusterMembership = zeros(A, N);

            obj.cluster_means = zeros(K, N);

            obj.cluster_errors = zeros(K, N);
            obj.cluster_sizes = zeros(K, N);
            obj.cluster_confidences = zeros(K, N);

            nSensors = expSettings.nSensors();
            obj.sensor_cluster_credit = zeros(nSensors,K);
            obj.sensor_cluster_weight = zeros(nSensors,K);
            obj.sensor_cluster_acc_credit = zeros(nSensors,K);
            
            obj.clusterData = expSettings.clusterData;
        end
        
        function cluster(obj, baseLogodds, basePost, labels, agentConfidence)
            
            import clustering.*
            
            A = obj.A; %no. agents
            K = obj.K; %no. clusters
            L = obj.L; %length of clustering window
            
            %Need to run these scripts first:
            %trainAndTestRun
            %cluster_agents_kmeans

            %cluster the agents
            %results
            rmpath(fullfile(pwd, '', 'netlab_3_3'));

            %first iteration
            cluster_weighter = dataReweighter(K, L);

            if strcmp(obj.clusterData, 'logodds')
                data = baseLogodds;
            elseif strcmp(obj.clusterData, 'post')
                data = basePost;
            else
                display('Defaulting to clustering on logodds: unrecognised setting');
                data = baseLogodds;
            end

            for n=1:obj.N

                display(sprintf('clustering iteration %i', n) );

                if n < L
                    start = 1;
                else
                    start = n - L + 1;
                end  

                if obj.weightClusterInput
                    cluster_input = cluster_weighter.reweight_data(data(:, 1:n) );
                else

                    cluster_input = data(:, start:n);
                end

                if n==1
                    clusters = kmeans( cluster_input, K, ...
                                    'EmptyAction', 'Singleton', 'Start', 'Uniform');
                else
                    %work out means for current sequence for previous clusters - use as
                    %seeds for current cluster
                    num_dim = L;
                    if n < L
                        num_dim = n;
                    end

                    seeds = zeros(K, num_dim );
                    %need to calculate mean across all weighted dimensions
                    for k=1:K
                        seeds(k, :) = sum(cluster_input(cluster_labels(clusters)==k, :));
                        seeds(k, :) =  seeds(k, :) ./ obj.cluster_sizes(k, n-1);
                    end

                    clusters = kmeans( cluster_input, K, ...
                        'EmptyAction', 'Singleton', ...
                        'Start', seeds );
                end  

                for k=1:K
                    current_clusterData = data(clusters==k, n);
                    obj.cluster_means(k, n) = sum(current_clusterData) / length(current_clusterData);
                    obj.cluster_sizes(k, n) = length(current_clusterData);
                end

                if obj.weightClusterInput
                    [cluster_labels] = labelClusters(clusters, ...
                        A, K, ...
                        cluster_weighter.reweight_data(obj.cluster_means(:,start:n)), obj.clusterData );        
                else
                    [cluster_labels] = labelClusters(clusters, ...
                        A, K, obj.cluster_means(:,start:n), obj.clusterData );
                end

                obj.cluster_means(cluster_labels, n) = obj.cluster_means(:, n);
                obj.cluster_sizes(cluster_labels, n) = obj.cluster_sizes(:, n);

                %need the mean of the posteriors not of the logodds to find the error
                for k=1:K
                    obj.cluster_errors(k, n) = ( sum(basePost(clusters==k,n)) ...
                        /obj.cluster_sizes(k,n) - labels(n,1) )^2;
                end

                for k=1:K
                    obj.cluster_confidences(k, n) = ...
                        sum(agentConfidence(cluster_labels(clusters)==k, n)) ...
                        / obj.cluster_sizes(k,n);
                end

                obj.agentClusterMembership(1:A, n) = cluster_labels(clusters);

            end
        end
    end
end