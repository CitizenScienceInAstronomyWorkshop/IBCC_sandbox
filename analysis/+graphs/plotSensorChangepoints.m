function plotSensorChangepoints( changes, simple )
%plots sensor changepoints        

    %data points where a sensor becomes informative
    obj.sensorPlusN = changes{1};
    %records which sensors become informative at each changepoint
    obj.sensorPlus = changes{2};
    obj.sensorMinusN = changes{3};
    obj.sensorMinus = changes{4};

    %draw the changepoint on the graph as well
    if simple
        obj.sensorPlus = ones(1, length(obj.sensorPlusN));
        obj.sensorMinus = ones(1, length(obj.sensorMinusN));
        plot(obj.sensorPlusN, obj.sensorPlus, '>', ...
            'MarkerSize', 10, 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0]);
        hold on

        plot(obj.sensorMinusN, obj.sensorMinus, ...
            '>', 'MarkerSize', 10, 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0]);
        hold on        
    else
        plot(obj.sensorPlusN, obj.sensorPlus, '^', ...
            'MarkerSize', 10, 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0]);
        hold on   
    
        plot(obj.sensorMinusN, obj.sensorMinus, ...
            'v', 'MarkerSize', 10, 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0]);
        hold on
    end
end

