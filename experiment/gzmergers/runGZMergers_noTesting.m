outputDirLabel = 'output_negLabs1';
settings.gz_mergers_notesting;

%idxsToKeep = 1:525560; subLabels = labels;
%pick the combination methods to use.
combMethods = {...
    combiners.MeanDecision.shortLabel,...  
    combiners.SimpleMajorityVoting.shortLabel, ...
...%     combiners.weighted.WeightedSum.shortLabel,...
...%    combiners.weighted.WeightedMajority.subShortLabel, ...
...%     combiners.bcc.DDIbccVb2.shortLabel,...
...%       combiners.bcc.DynIbccVb.shortLabelInitAll,...
      combiners.bcc.IbccVb.shortLabel,...       
...%     combiners.bcc.IbccVb.shortLabelSeq, ...      
...%     combiners.bcc.Ibcc.shortLabel, ...
    };

sepPosLabFile = [expSet.getDataDir() '/positiveLabels.csv'];
sepNegLabFile = [expSet.getDataDir() '/negativeLabelsEasy.csv'];
rawColumnMap = [0 3 1 4 0 0];
columnString = '%u64 %u64 %u64 %u64';

clear idxsToKeep pointsToDrop baseOutputs rawData c_class

runExperiment
resultsAllFolds1 = resultsAllFolds;
testIdxs1 = testIdxsAllFolds;

auc1 = auc;
aucInt1 = aucInt;
lam1 = lam;

outputDirLabel = 'output_negLabs2';
settings.gz_mergers_notesting;

display('running with more labels - clearing everything first!!!!!!!!!!!!!');
clear idxsToKeep pointsToDrop baseOutputs rawData c_class
sepNegLabFile = [expSet.getDataDir() '/negativeLabelsTogether.csv'];
rawColumnMap = [0 3 1 4 0 0];
columnString = '%u64 %u64 %u64 %u64';
runExperiment

auc2 = auc;
aucInt2 = aucInt;
lam2 = lam;