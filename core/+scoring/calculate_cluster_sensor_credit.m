function [ sensor_acc_credit ] = calculate_sensor_credit( agents, nAgents, nSensors, nSamples, classifications, results, clusters, k )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    sensor_acc_credit = zeros(nSensors, nSamples);
        
    for n=1:nSamples
        for a=1:nAgents
            agent = agents{a};
            
            if clusters(a, 1) ~= k
                continue
            end
            
            classification = classifications(n, 1);
            if classification ~= -1
                logit = results(a, n);
            
                agent_credit = calculate_agent_sensor_credit(agent, classification, logit);
            
                sensor_acc_credit(:, n) = sensor_acc_credit(:, n) + agent_credit;
            end
        end
        
        if n>1
            sensor_acc_credit(:, n) =  sensor_acc_credit(:, n) + sensor_acc_credit(:, n-1) .* 0.9;
        end
    end
end

