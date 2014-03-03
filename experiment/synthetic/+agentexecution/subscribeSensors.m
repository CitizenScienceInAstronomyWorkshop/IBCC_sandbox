function [ subscriptions ] = subscribeSensors( nAgents, nSensors, expSettings)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    maxSensors = expSettings.maxSensors;
    minSensors = expSettings.minSensors;

    subscriptions = {nAgents};
    
    if expSettings.nInfSensors == 1 && expSettings.nNoninfSensors == 1 && nAgents == 5
        subscriptions{1} = [1 0];
        subscriptions{2} = [0 1];
        subscriptions{3} = [1 1];
        subscriptions{4} = [0 1];
        subscriptions{5} = [1 0];
        return;
    end
    
    if expSettings.nInfSensors==expSettings.nInformedAgents && ...
            expSettings.nNoninfSensors==expSettings.nAgents-expSettings.nInformedAgents
        for a=1:expSettings.nAgents
            subscriptions{a} = zeros(1, expSettings.nInfSensors+expSettings.nNoninfSensors);
            subscriptions{a}(a) = 1;
        end
        return;
    end
    

    all_sensors_list = {nSensors};
    for i=1:nSensors
        all_sensors_list{i} = i;
    end    
    
    
    
    for i=1:nAgents
        num_agent_sensors = (maxSensors - minSensors + 1) * rand(1) + minSensors;
        num_agent_sensors = floor(num_agent_sensors);
        sensors = ones(1, nSensors); %list of all sensors
        
        if i > expSettings.nInformedAgents && ...
                    expSettings.nInformedAgents > 0 
                
                sensors(1, 1:expSettings.nInfSensors) = 0;
                
        end
        while sum(sensors) > num_agent_sensors
            
           
            s_to_remove = ceil(rand(1) * numel(sensors));
            if i > expSettings.nInformedAgents || ...
                    sum(sensors(1:expSettings.nInfSensors)) > 1 || ...
                    s_to_remove > expSettings.nInfSensors
                sensors(1, s_to_remove) = 0;
            end
        end
        subscriptions{i} = sensors;
    end
end

