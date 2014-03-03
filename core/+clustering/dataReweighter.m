classdef dataReweighter
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        weights %for emphasizing more recent data points
        k %number of clusters
        seqLength %number of dimensions or (values from one item) by which to cluster
    end
    
    methods
        function obj= dataReweighter(k, seqLength)
            obj.k = k;
            obj.seqLength = seqLength;
            
            %initialize weights
            obj.weights = zeros(seqLength, 1);
            for d=1:seqLength
                obj.weights(d, 1) = (1.0+d) / (1.0+seqLength);
            end
        end
        
        function reweighted = reweight_data(obj, phi)
            
            data_size = size(phi);
            
            %number of dimensions available in the data
            seqLength_available = data_size(2);
            
            %number of dimensions we are going to use (we used the last D
            %dimensions)
            seqLength = obj.seqLength;
            
            if seqLength > seqLength_available
                seqLength = seqLength_available;
            end
            
            reweighted = phi(:, seqLength_available-seqLength+1:seqLength_available);
            reweighted = reweighted .* ...
                transpose(obj.weights(obj.seqLength-seqLength+1:obj.seqLength, 1) * ones(1, data_size(1)));
                        
        end
        
        function means = init_means(phi)
            %pick k different values as from the data set as cluster
            %centres. If there are < k different values in the dataset,
            %lower the value of k accordingly
            
            means = zeros(obj.k, 1);
            inited = 0;
            seed = 1;
            
            num_seeds = size(phi);
            num_seeds = num_seeds(2)            
            
            while inited < obj.k && seed <= num_seeds
                mean = phi(seed);
                means(inited) = mean;
                
                inited = inited + 1;
                seed = seed + 1;
                
                %check we don't have this value already
                for i=1:inited
                    if means(i) == mean
                        %put the number back so next iteration is the same
                        inited = inited - 1; 
                        break
                    end
                end
                
                means
            end
            
            
        end
    end
    
end

