function [ colourSet ] = createColourSet( nSensors )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    colourSet = zeros(nSensors, 3);


    if nSensors <= 7
        colourSet = [0 0 1; 1 0 0; 0 1 0; 0 0 0; 1 1 0; 0 1 1; 1 0 1];
        return;
    end

    %avoid white out
    whiteLimit = 7/8;

    base = (nSensors / whiteLimit) ^ (1/3);

    base = ceil(base);

    for s=1:nSensors

    %     colour = s / nSensors;
    %     
    %     
    %     colourSet(s, 1) = s * 4;
    %     colourSet(s, 1) = mod(colourSet(s, 1), nSensors); 
    %     colourSet(s, 1) = colourSet(s, 1) / nSensors;
    %     
    %     colourSet(s, 3) = colour * 2;
    %     if(colourSet(s, 3) > 1)
    %         colourSet(s, 3) = colourSet(s, 3) -1 ;
    %     end
    %     
    %     colourSet(s, 2) = 1 - colour;
    %     
    %     if sum(colourSet(s, :)) > 2
    %         
    %     end

        key = s-1;

        colourSet(s, 1) = mod(key, base);
        colourSet(s, 2) = floor(mod(key, base^2)/base);
        colourSet(s, 3) = floor(key/base^2);

        colourSet(s, :) = colourSet(s, :) / (base-1);
    end