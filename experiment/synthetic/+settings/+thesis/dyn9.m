function [expSet, synthSet, bccSet, rerun, range, varName] = dyn9(i)

    varName = 'Slow changes over 200 data points, more evenly spaced';%'Mean Individual AUC';
    rerun = true; display('warning: likely to squash any saved datasets');
    settings.thesis.dyn0;
    expSet.setDirNames([rootForAllExpts '/exp9_moresmoothing_rnd100_2/']);

    range = 0;%[0 0.2];
    synthSet.switchDuration = 100;
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
    synthSet.defNoiseChanges = sort(randi(800, synthSet.nDatasets, 15)')' + 100;%[100 200 280 340 400 460 540 620 700 800];
end