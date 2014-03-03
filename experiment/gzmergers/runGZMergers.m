
settings.gz_mergers1;

%idxsToKeep = 1:525560; subLabels = labels;
%pick the combination methods to use.
combMethods = {...
    combiners.MeanDecision.shortLabel,...  
    combiners.SimpleMajorityVoting.shortLabel, ...
...%     combiners.weighted.WeightedSum.shortLabel,...
...%    combiners.weighted.WeightedMajority.subShortLabel, ...
...%     combiners.bcc.DDIbccVb2.shortLabel,...
...%       combiners.bcc.DynIbccVb.shortLabelInitAll,...
      combiners.bcc.IbccDiff.shortLabel,...
...%      combiners.bcc.IbccPooledVb.shortLabel,...       
      combiners.bcc.IbccVb.shortLabel,...       
...%     combiners.bcc.IbccVb.shortLabelSeq, ...      
...%     combiners.bcc.Ibcc.shortLabel, ...
    };

sepPosLabFile = [expSet.getDataDir() '/positiveLabels.csv'];
sepNegLabFile = [expSet.getDataDir() '/negativeLabelsTogether.csv'];
rawColumnMap = [0 3 1 4 0 0];
columnString = '%u64 %u64 %u64 %u64';

%clear idxsToKeep pointsToDrop baseOutputs rawData c_class

runExperiment

auc1 = auc;
aucInt1 = aucInt;
lam1 = lam;
labelsBoth = subLabels;
% 
% display('running with different labels - clearing everything first!!!!!!!!!!!!!');
% clear idxsToKeep pointsToDrop baseOutputs rawData c_class
% sepNegLabFile = [expSet.getDataDir() '/negativeLabelsEasy.csv'];
% rawColumnMap = [0 3 1 4 0 0];
% columnString = '%u64 %u64 %u64 %u64';
% runExperiment
% 
% auc2 = auc;
% aucInt2 = aucInt;
% lam2 = lam;
% 
% aucInt2_allLabs = aucIntAlt;


%  display('running with more labels - clearing everything first!!!!!!!!!!!!!');
%  clear idxsToKeep pointsToDrop baseOutputs rawData c_class
%  sepNegLabFile = [expSet.getDataDir() '/negativeLabelsHard.csv'];
%  rawColumnMap = [0 3 1 4 0 0];
%  columnString = '%u64 %u64 %u64 %u64';
%  runExperiment
%  
%  auc3 = auc;
%  aucInt3 = aucInt;
%  lam3 = lam;