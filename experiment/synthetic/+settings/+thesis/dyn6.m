function [expSet, synthSet, bccSet, rerun, range, varName] = dyn6(i)

    varName = 'Lots of Evenly-spaced, sudden changes';%'Mean Individual AUC';
    rerun = true; display('warning: likely to squash any saved datasets');
    settings.thesis.dyn0;
    expSet.setDirNames([rootForAllExpts '/exp6/']);

    range = 0;%[0 0.2];
    
    bccSet.fixedNIt = 0;
    
%     synthSet.means = synthSet.means .* range(i);
%     synthSet.deviations = synthSet.deviations .* 0.1;%range(i).^0.5;
    synthSet.pMode = [0.5 0 0.5 0 0.5 0 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 ];%0.5 0.5 0.5 0.5 0.5];% 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5];%0.5 0.5 0 0.5];%[0.5 range(i) 0.5 range(i) 0.5 range(i) 0.5 range(i) 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5];
%     synthSet.deviations = [synthSet.deviations repmat(synthSet.deviations(:,end), 1, synthSet.nAgents-5)];
%     synthSet.means = [synthSet.means repmat(synthSet.means(:,end), 1, synthSet.nAgents-5)];
    synthSet.nInfSensors = synthSet.nAgents;
    synthSet.nNoninfSensors = 0;
    synthSet.maxSensors = 1;
    synthSet.minSensors = 1;
    synthSet.nInformedAgents = synthSet.nAgents;
    synthSet.class1ErrorOnly = 0;
    synthSet.defNoiseChanges = [100 200 280 340 400 460 540 600 660 720 800 900];%[140 280 420 560 700 840];
end