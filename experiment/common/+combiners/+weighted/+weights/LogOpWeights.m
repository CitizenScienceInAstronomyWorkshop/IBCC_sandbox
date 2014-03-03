classdef LogOpWeights < combiners.weighted.weights.Weights
    %PRECISIONWEIGHTS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        agents        
        bias
    end
    
    methods
        
        function obj = LogOpWeights(agents)
            obj.agents = agents; 
        end
        
        function logodds = getLogodds(obj, post)
            %base classifiers may produce 0s and 1s that are wrong. These
            %are probably mis-calibrations, so let's not exclude them;
            %recalibrate.
            negPost = post<0.5;
            post = abs(post - 0.5).^1.01 + 0.5;
            post(negPost) = 1 - post(negPost);            
            
            logodds = log(post) - log(1-post);
        end
        
        function W = calculateWeights(obj, post, ~, T)           

            logodds = obj.getLogodds(post);
            
            mu = obj.calculateMu(logodds,  T);
            sigmaSq = obj.calculateSigSq(logodds, mu);
                 
            
            W = 2 .* mu ./ (1+sigmaSq);
        end
        
        function mu =  calculateMu(obj, logodds,  T)
            %calculate mean for each class
            
            nPos = (sum(T==1, 2));
            nNeg = (sum(T==0, 2));
                        
            if nPos==0 || nNeg==0
                obj.bias = zeros(size(logodds));
                mu = zeros(size(logodds,1),1);
            else
                posMean = (sum(logodds(:, T==1), 2) ) ./ nPos;
                negMean = -(sum(logodds(:, T==0), 2)) ./ nNeg;                
                obj.bias = (posMean - negMean) ./2;
                obj.bias = repmat(obj.bias, 1, size(logodds,2));
                mu = (posMean + negMean) ./ 2;
            end
        end
        
        function var = calculateSigSq(obj, logodds, mu)            
            logodds = logodds - obj.bias;
            diffs = (abs(logodds) - repmat(mu, 1, size(logodds,2))).^2;
            var = sum(diffs, 2) ./ size(diffs,2);
            
        end
             
    end
    
end

