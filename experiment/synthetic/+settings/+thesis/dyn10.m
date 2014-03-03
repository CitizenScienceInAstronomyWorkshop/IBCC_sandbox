function [expSet, synthSet, bccSet, rerun, range, varName] = dyn10(i)

    varName = 'Slow changes over 200 data points, randomly spaced';%'Mean Individual AUC';
    rerun = true; display('warning: likely to squash any saved datasets');
    settings.thesis.dyn0;
    expSet.setDirNames([rootForAllExpts '/exp10/']);

    range = 0;%[0 0.2];
    synthSet.switchDuration = 200;
    bccSet.fixedNIt = 0;

%     synthSet.means = synthSet.means .* range(i);
%     synthSet.deviations = synthSet.deviations .* 0.1;%range(i).^0.5;
    synthSet.pMode = [0.5 0 0.5 0 0.5 0 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 ];
    synthSet.class1ErrorOnly = 0;
    synthSet.nInfSensors = synthSet.nAgents;
    synthSet.nNoninfSensors = 0;
    synthSet.maxSensors = 1;
    synthSet.minSensors = 1;
    synthSet.nInformedAgents = synthSet.nAgents;
    synthSet.defNoiseChanges = randi(1000, synthSet.nDatasets, 10);
end