function [ informedness ] = plotInformedness( agent, nInfSensors, nSamples, ...
    changes, min, max)
% Plot a graph of how many informed sensors an agent has
% Informedness is a table: 1st column is n-value (datapoint) of a change or 
% last data point at a given informedness value.
% 2nd column is the new value of informedness

%     %informedness lists x values of changepoints; y values of 0 for when
%     %agent becomes uninformed, y values for no. informative sensors it has
%     informedness = [];
%     informedness(1, 1) = 1; %1st data point
%     informedness(1, 2) = 0; %assume uninformed to start    
%     
%     for s=1:nInfSensors
%         if agent.s(1, s) == 1
%             informedness(1, 2) = informedness(1, 2) + 1; 
%             % +1 for no. informative sensors
%         end
%     end
%     
%     
%     num_changepoints = 1;
%     %skip the inf sensors at start as we've already done these
%     for change_idx=nInfSensors+1:length(changepoints_x)
%         s_minus = changepoints_minus(1, change_idx);
%         if agent.s(1, s_minus) == 1
%             %less informed!
%             
%             %last point at the old value
%             num_changepoints = num_changepoints + 1;           
%             
%             informedness(num_changepoints, 1) = ...
%                 changepoints_x(1, change_idx) - 1;
%             
%             informedness(num_changepoints, 2) = ...
%                 informedness(num_changepoints -1, 2);
%             
%             %new value
%             num_changepoints = num_changepoints + 1;
%             
%             informedness(num_changepoints, 1) = ...
%                 changepoints_x(1, change_idx);
%             
%             informedness(num_changepoints, 2) = ...
%                 informedness(num_changepoints -1, 2) - 1;
%             
%         end
%         
%         s_plus = changepoints_plus(1, change_idx);
%         if agent.s(1, s_plus) == 1
%             %last point at old value
%             
%             num_changepoints = num_changepoints + 1;
%       
%             informedness(num_changepoints, 1) = ...
%                 changepoints_x(1, change_idx) - 1;
%             
%             informedness(num_changepoints, 2) = ...
%                 informedness(num_changepoints -1, 2);
%             
%             %new value
%             num_changepoints = num_changepoints + 1;
%       
%             informedness(num_changepoints, 1) = ...
%                 changepoints_x(1, change_idx);
%             
%             informedness(num_changepoints, 2) = ...
%                 informedness(num_changepoints -1, 2) + 1;
%         end
%     end
%     num_changepoints = num_changepoints + 1;
%     informedness(num_changepoints, 1) = nSamples;
%     informedness(num_changepoints, 2) = informedness(num_changepoints-1, 2); 
%     
    %NEW!
    
    %1. Get the changepoints that apply to this agent
    changesPlusN = changes{1};
    changesPlusS = changes{2};
    
    relevantPlus = ismember(changesPlusS, find(agent.s));
    changesPlusN = changesPlusN(relevantPlus);
    
    changesMinusN = changes{3};
    changesMinusS = changes{4};
    
    relevantMinus = ismember(changesMinusS, find(agent.s));
    changesMinusN = changesMinusN(relevantMinus);
    
    %2. Sum up the pluses and minuses
    
    %get the points immediately after a change to the sensors
    changesCombined = union(changesMinusN, changesPlusN);
    if isempty(changesCombined) || changesCombined(1) ~= 1
        changesCombined = [1 changesCombined];
    end
    %add in the points immediately before each change
    changesCombined = [changesCombined (changesCombined-1)];
    
    if changesCombined(length(changesCombined)) ~= nSamples
        changesCombined = [changesCombined nSamples];
    end
    
    changesCombined = sort(changesCombined);
    informedness = zeros(1, nSamples);
    
    for i=1:length(changesPlusN)
        informedness(changesPlusN(i):nSamples) = informedness(changesPlusN(i):nSamples) + 1;
    end

    for i=1:length(changesMinusN)
        informedness(changesMinusN(i):nSamples) = informedness(changesMinusN(i):nSamples) - 1;
    end  
    
    %get the sample indices for when a change to the sensors occurs. The 
    %informedness value will not necesarily change but it is useful to mark
    %the changepoint on the graph.
    informedness = [informedness(1) informedness];
    informedness = informedness(changesCombined+1);
    
    hold all
    plot( changesCombined, informedness/nInfSensors * (max-min) + min, ... 
        'x-', 'Color', [0.5 0.5 0], 'LineWidth', 1);
        
end

