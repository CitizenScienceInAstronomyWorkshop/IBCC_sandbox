colourSet = createColourSet(nSensors);

%plot the sensor credit and sensor weight values
figure(5);

legendStrings = {};
for s=1:nSensors
    plot( 1:nSamples, sensor_credit(s, :), 'Color', colourSet(s, :));
    title('Sensor credit')
    hold all
    
    sensor_str = sprintf('Sensor %d', s);
    legendStrings{s} = sensor_str;
end

legend( char(legendStrings) );

hold off

%plot the sensor credit and sensor weight values
figure(6);

for s=1:nSensors
    %if(colour > 7)
    %    colourSet(s, 3) = (colour - 7);
    %end
    
    plot( 1:nSamples, sensor_weight(s, :), 'Color', colourSet(s,:) );
    title('Sensor weight credit')
    hold all
end

legend( char(legendStrings) );

hold off

%plot the sensor credit and sensor weight values
figure(7);

for s=1:nSensors

    plot( 1:nSamples, sensor_acc_credit(s, :), 'Color', colourSet(s,:) );
    title('Accumulating Sensor Credit')
    hold all
end

legend( char(legendStrings) );

hold off