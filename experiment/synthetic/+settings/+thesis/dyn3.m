function [expSet, synthSet, bccSet, rerun, range, varName] = dyn3(i)

    varName = 'Sensor Error Rate';%'Mean Individual AUC';
    rerun = true; display('warning: likely to squash any saved datasets');
    settings.thesis.dyn0;
    expSet.setDirNames([rootForAllExpts '/exp3_moresmoothing/']);

    range = [0 0.1];
    synthSet.switchDuration = 1;
    synthSet.defNoiseChanges = [];
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

end