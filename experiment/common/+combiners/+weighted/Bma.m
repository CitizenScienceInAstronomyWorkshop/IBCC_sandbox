classdef  Bma < combiners.weighted.WeightedSum
    %Bma Bayesian Model Averaging. Performs soft selection rather than
    %combination but it will be interesting to see how it filters out
    %uninformed agents
    
    properties        
        pk %prior probability of classifiers
        nPriors = 5 %amount of made up data that biases the weights towards the prior probabilities
        
        W = []
        members
    end
    properties (Constant)
        subSubShortLabel = 'bma';
    end
    
    methods
        function obj = Bma(nAgents, K, targets, agents, windowSize)
            obj@combiners.weighted.WeightedSum(nAgents, K, targets, agents, windowSize);
            obj.pk = 1/nAgents;
            obj.label = 'Bayesian Model Averaging';
        end     
        
        function combinedPost = combineCluster(obj, post, T, members, k)
            if ~exist('k','var')
                k = 1;
            end
            combinedPost = obj.combineCluster@combiners.weighted.WeightedSum(post, T, members, k);
            %obj.plotWeights(obj.W, k);
        end

        function scores = getScores(obj, post)
            if obj.discreetInputs
                scores = round(post);
            else
                scores = post;
            end
        end
        
        function createWeightCalculator(obj, members, ~)
            obj.W = [];
            obj.members = members;
            obj.weightCalculator = obj;
        end
        
        function W = calculateWeights(obj, ~, post, T, ~, n)
                
            %ignoring any data points we don't have labels for
            post = post(:, T(1,:)~=-1);
            T = T(:, T(1,:)~=-1);
            
            %add some data to the start that balances the weights initially
            start = n - obj.windowSize;
            if start < 1
                start = 1;
            end
            
            if n-1 > size(T,2)
                iEnd = size(T,2);
            else
                iEnd = n-1;
            end
            
            post = [ones(obj.nAgents, obj.nPriors) post(:, start:iEnd)];
            %base classifiers may produce 0s and 1s that are wrong. These
            %are probably mis-calibrations, so let's not exclude them;
            %recalibrate.
            negPost = post<0.5;
            post = abs(post - 0.5).^1.01 + 0.5;
            post(negPost) = 1 - post(negPost);
                        
            T = [ones(1, obj.nPriors) T(start:iEnd)];
            
            %k - base classifier
            %D - data 
            %T - part of the data; known targets!!! This probably fails if they're not supplied
            %X - part of the data; posteriors from the agents
            %yi - prediction for data point i by k
            % ti - target label for i
            
            %p(k|D) = p(D|k)p(k) / p(D)
                        
            %find model evidence p(D|k) = p(X, T|k)
            % = p(T|k,X) * p(X|k)
            % = prod_i( yi^ti * (1-yi)^(1-ti) ) * p(X,k)/p(k) %%% first term is probability of the target labels according to k; second term is 1 as we observe it?
            
            T = ones(obj.nAgents, 1) * T;
            
                           
            pti_given_kx = post.^T .* (1-post).^(1-T); %first term above
            %don't count pCorrect for new members. Instead assign them an
            %average pCorrect of all other members.
            
            if ~isempty(obj.members)
                M = [ones(obj.nAgents, obj.nPriors) obj.members(:, start:iEnd)];            
                pti_given_kx = pti_given_kx .* M; %set irrelevant values to zero
                A = sum(M,1);
            else
                A = obj.nAgents;
                M = 0;
            end
            
            %for any non-members of the cluster, insert the average as an estimate
            avgPti_given_kx = sum( pti_given_kx, 1) ./ A; %average probability across the  agents
            estimates = ones(obj.nAgents, 1) * avgPti_given_kx;
            estimates = estimates .* (1-M);
            pti_given_kx = pti_given_kx + estimates;
            
            %completes the first term
            pTgivenkX = prod(pti_given_kx, 2);
            
            %priors
            if  ~isempty(obj.members)
                pk = obj.pk .* obj.members(:, n);
            else
                pk = obj.pk;
            end
            
            %p(D) = p(X, T) ->we don't need this to calculate the relative
            %probabilities of each classifier as the p(k|D)s sum up to one.
            % p(D) = sum_k(p(D|k)p(k))
            pTgivenX = sum(pTgivenkX.*pk, 1);
            
            %W = p(k|D)
            pTkgivenX = pTgivenkX .* pk;
            
            pkgivenTX = pTkgivenX ./ pTgivenX;
            Wn = pkgivenTX;
            obj.W = [obj.W Wn];
            
            W = obj.W;       
        end
       
    end
end

