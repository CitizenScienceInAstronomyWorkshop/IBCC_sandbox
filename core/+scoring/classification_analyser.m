classdef classification_analyser < handle
    
    properties
        
        data_similarity
        results
        classifications
        data_sim_total
        
    end
    
    methods
        
        function obj = classification_analyser(results, classifications, test_data)
            
            obj.calculate_data_similarity(test_data);
            
            obj.results = results;
            obj.classifications = classifications;

        end
        
        function calculate_data_similarity(obj, test_data)
           
            N = size(test_data, 1);
            
            obj.data_sim_total = zeros(1, N);
            obj.data_similarity = zeros(N, N);
            
            for n=1:N
                for m=1:n-1
                    obj.data_similarity(m, n) = exp(- sum(abs( (test_data(n, 1) - test_data(m, 1)) )) );
                    obj.data_sim_total(1,n) = obj.data_sim_total(1, n) + obj.data_similarity(m, n);
                end
            end
        end
        
        function [agentConfidence] = agentConfidence(obj, A, n)
            agentConfidence = zeros(length(A), 1);

            for m=1:n-1
                
                result_similarity = 1 - ...
                    abs( (ones(length(A),1)*obj.classifications(m, 1)) ...
                    - obj.results(A, m) );
                
                agentConfidence = agentConfidence + ...
                    ( ones(length(A),1)*obj.data_similarity(m, n)) .* result_similarity;
            end

            agentConfidence = (agentConfidence ./ obj.data_sim_total(1,n));
        end
    end
end

