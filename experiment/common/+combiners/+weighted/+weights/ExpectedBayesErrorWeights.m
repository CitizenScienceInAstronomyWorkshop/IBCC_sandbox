classdef ExpectedBayesErrorWeights < combiners.weighted.weights.Weights
    %EXPECTEDBAYESERROR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        nAgents
    end
    
    methods
        function obj = ExpectedBayesErrorWeights(nAgents)
            obj.nAgents = nAgents;
        end
        
        function W = calculateWeights(obj, scores, T)
            nSamples = length(logodds);            
            beta = zeros(obj.nAgents, nSamples);            
            for a=1:obj.nAgents

                %prior expectation of what the difference between the
                %output of an agent and the target will be
                prior_diff = 0.3;
                %how much weight do we give to the prior? Treat it as a
                %number of observations prior_weight, where the difference is
                %prior_diff
                prior_weight = 5; 

                diffs = abs(scores(a, :) - T);
                for n=nSamples:2

                    start = n-1-obj.windowSize;
                    if start < 1
                        start = 1;
                    end

                    currentWindowSize = n - start - 1;

                    diffs(n) = (sum(diffs(start:n-1)) + prior_diff*prior_weight) ...
                        / (currentWindowSize+prior_weight);
                end
                diffs(1) = prior_diff;
                %if diff < 0.5 then we expect this agent to be helpful,
                %so it is weighted positively. If diff==0.5 then we
                %ignore the agent. If diff > 0.5 then we expect the
                %opposite of what the agent thinks so assign a negative
                %weight.
                beta(a, :) = 0.5 - diffs;
            end     
            
            W = beta ./ (ones(obj.nAgents, 1) * sum(beta, 1));
        end
            
    end
    
end

