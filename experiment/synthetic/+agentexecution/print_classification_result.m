function print_classification_result(agent, logit, nInfSensors)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    inf_str = '';
    for s=1:nInfSensors
        if agent.s{s} <= nInfSensors
            inf_str = [inf_str '%d, '];
            inf_str = sprintf(inf_str, agent.s{s});
        else
            break
        end

    end

    fprintf('%d \t inf sensors=[%s] \t logit: %f\n', agent.label, inf_str, logit);

end

