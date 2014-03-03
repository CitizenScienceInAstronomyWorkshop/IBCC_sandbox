close all
clear progressMonitor
settings.thesis.gz_sn_subs; %works pretty well 92% versus 75% - just talk about balancing the training data; nothing about filtering agents
% settings.thesis.gz_sn_subs2; %original setup from posters
% settings.thesis.gz_sn_subs3; %more fiddling. Try to close in on poster-level performance with a more plausible subsampling technique

%idxsToKeep = 1:525560; subLabels = labels;
%pick the combination methods to use.
combMethods = {...
    combiners.MeanDecision.shortLabel,...  
    combiners.SimpleMajorityVoting.shortLabel, ...
    combiners.weighted.WeightedMajority.subShortLabel, ...
 ...%     combiners.bcc.DynIbccVb.shortLabelInitAll,...
 ...%combiners.bcc.IbccDiff.shortLabel,...         
    combiners.bcc.IbccVb.shortLabel,...       
    combiners.bcc.IbccSampling.shortLabel, ...
    combiners.weighted.WeightedSum.shortLabel,...
    };
legalScores = [-1 1 3];
%not sure if the interpretation of screen labels is helpful.

%Test how well we do when empty screen labels are interpreted as negatives
%in the test labels but are not provided for training. I.e. how well we do
%in the real system over these objects. At the moment when we ignore the
%empty screen labels we do poorly, but it's possible we're still doing well
%on the remaining unlabelled objects that we haven't checked.

nRounds = 4;
nNegsPerRound = floor(1091/nRounds);
totalPerRound = (1091+330)/4;

resultsAllRounds = [];
labelsAllRounds = [];
avgAucAllRounds = 0;
brierAllRounds = 0;

for round=1:nRounds
    
    expSet.startNegAssets = 1 + (round-1)*nNegsPerRound;
    
    runExperiment
    
    avgAucAllRounds = avgAucAllRounds + avgAuc;
    
    resultsAllRounds = [resultsAllRounds resultsAllFolds(:, testIdxsAllFolds)];
    labelsAllRounds = [labelsAllRounds labelsAllFolds(testIdxsAllFolds)];
    brierAllRounds = brierAllRounds + brier;
    pointsToDrop = [];
    clear idxsToKeep
end

close all

%DRAW A LOVELY ROC CURVE --------------------------------------------------
auc = graphs.ClassifierPerformanceGraph.drawRoc(resultsAllRounds, labelsAllRounds-1, combMethods, false, false, true)

avgAucAllRounds ./ nRounds
brierAllRounds = brierAllRounds ./ nRounds