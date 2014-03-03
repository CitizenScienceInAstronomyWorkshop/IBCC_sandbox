function [expSet, synthSet, bccSet, rerun, range, varName] = dyn4(i)

    varName = 'Slow changes over 200 data points';%'Mean Individual AUC';
    rerun = true; display('warning: likely to squash any saved datasets');
    settings.thesis.dyn0;
    expSet.setDirNames([rootForAllExpts '/exp4/']);

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
    synthSet.defNoiseChanges = [100 200 250 300 350 400 650 700 750 800 ];
end