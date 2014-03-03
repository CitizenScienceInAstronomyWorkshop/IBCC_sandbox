function [expSet, synthSet, bccSet, rerun, range, varName] = exp2(i)

    varName = 'Class 1 Error Rate for Agents 1 to 3';%added equally to all sensors: strength of over-exaggeration and under-exaggeration of confidence';
    rerun = true; display('warning: likely to squash any saved datasets');
    settings.thesis.exp0;
    expSet.setDirNames([rootForAllExpts '/exp2_mvredo/']);
    
    %data generation
%     synthSet.means = [synthSet.means ];%[1;-1]];
%     synthSet.deviations = synthSet.deviations .* 0.1;% [2; -2]];
%     synthSet.nAgents = synthSet.nAgents + 1;
%     synthSet.nInformedAgents = synthSet.nAgents;
%     synthSet.nInfSensors = synthSet.nInformedAgents;
%     synthSet.cal = [synthSet.cal 1];
%     synthSet.bias = [synthSet.bias 0];
    range = [0.1 0.3 0.5 0.7 0.9]; %working for case where we lose
%     information as the value goes up.
    
    synthSet.class1ErrorOnly = 1;
    synthSet.pMode = zeros(1, synthSet.nAgents);
%     synthSet.pMode(1:3) = synthSet.pMode(1:3) + range(i);
    synthSet.pMode = synthSet.pMode + range(i);

%     synthSet.deviations(1,:) = synthSet.deviations(1,:) + range(i);
end