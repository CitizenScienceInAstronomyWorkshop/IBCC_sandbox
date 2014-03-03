function [ agent_credit ] = calculate_agent_sensor_credit( agent, classification, logit )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


    agent_credit = (0.5 - abs(classification-logit)) ./ sum(agent.s);

    if(agent_credit < 0)
        agent_credit = agent_credit .* transpose(agent.s)* 3; %punish errors
    else
        agent_credit = agent_credit .* transpose(agent.s);
    end

end

