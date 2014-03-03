function [ D T ] = generateTwoClassData( expSettings )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    p_c1 = expSettings.p_c1;
    means = expSettings.means;
    deviations = expSettings.deviations;
    nInfSensors = expSettings.nInfSensors;
    nSensors = expSettings.nSensors();
    nTrainingSamples = expSettings.nTrainingSamples*expSettings.nAgents;
    
    %pMode = expSettings.pMode;

    D = randn(nTrainingSamples, nSensors);
    T = zeros([ nTrainingSamples, 1]); %target values

    for n=1:nTrainingSamples
        [ row, t ] = datageneration.generateDataPoint(...
            nSensors, nInfSensors, means, deviations, p_c1, 0);

        T(n) = t;
        D(n,:) = row;
    end
end

