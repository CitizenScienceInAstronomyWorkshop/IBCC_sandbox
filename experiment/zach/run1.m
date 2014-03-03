
settings.zach1;

combMethods = {...
    combiners.MeanDecision.shortLabel,...  
    combiners.SimpleMajorityVoting.shortLabel, ...
    combiners.bcc.IbccVb.shortLabel,...       
    };

rawColumnMap = [0 1 2 3 0 4];
columnString = '%u64 %u64 %d %d';
legalScores = [-1, 1];

% clear idxsToKeep pointsToDrop baseOutputs rawData c_class

runExperiment
resultsAllFolds1 = resultsAllFolds;
testIdxs1 = testIdxsAllFolds;

auc1 = auc;
aucInt1 = aucInt;
lam1 = lam;