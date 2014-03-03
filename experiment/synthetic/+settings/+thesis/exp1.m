function [expSet, synthSet, bccSet, rerun, range, varName] = exp1(i)

    varName = 'Sensor Error Rate';%'Mean Individual AUC';
    rerun = true; display('warning: likely to squash any saved datasets');
    settings.thesis.exp0;
    expSet.setDirNames([rootForAllExpts '/exp1_mvredo/']);

%     range = [16 4 1 0.5.^2 0.25.^2];%[0 0.5 1 2 4];
    range = [0 0.1 0.2 0.3 0.4 0.45 0.48 0.5];%
    
    bccSet.fixedNIt = 200;

%     synthSet.means = synthSet.means .* range(i);
%     synthSet.deviations = synthSet.deviations .* 0.1;%range(i).^0.5;
    synthSet.pMode = zeros(1,synthSet.nAgents) + range(i);
    synthSet.class1ErrorOnly = 0;
end