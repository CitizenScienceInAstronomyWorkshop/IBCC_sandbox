function [expSet, synthSet, bccSet, rerun, range, varName] = dyn5(i)

    varName = 'Evenly-spaced, sudden changes';%'Mean Individual AUC';
    rerun = true; display('warning: likely to squash any saved datasets');
    settings.thesis.dyn0;
    expSet.setDirNames([rootForAllExpts '/exp5_moresmoothing_rnd100_2/']);

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
    synthSet.defNoiseChanges = sort(randi(800, synthSet.nDatasets, 15)')' + 100; %avoid changes at very beginning and end as they wouldn't affect results much%[100 200 300 400 460 540 600 700 800 900];%[140 280 420 560 700 840];
end