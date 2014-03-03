%Run this code from this directory
%clear
import settings.*

%load('/homes/49/edwin/results/agent_clustering_meam-magnitude_ordered_relabelling/data/5.mat');

bccSettingsSimple1

%avoid clashes with multi-run tests
expSettings.topSaveDir = sprintf('%s/singleRun', expSettings.topSaveDir);

%cluster agents over a data sequence of 20. use posteriors as cluster input
%data.
% expSettings.clusterData = 'post';
% expSettings.weightClusterInput = false;
% expSettings.topSaveDir = '/homes/49/edwin/results/combination/inf5_noninf11_agents5dlr';
% expSettings.expLabel = 'dlrCombiners';
%reuseSubs = true;

%pick the combination methods to use.
combMethods = {...
    combiners.MeanDecision.shortLabel, ...
...%     combiners.bcc.IbccVb.shortLabelSeq,...
...%    combiners.bcc.DynIbccVb.shortLabelInitAll, ...
    combiners.bcc.DynIbccVb.shortLabel, ...
    combiners.bcc.IbccVb.shortLabel, ...
%     combiners.weighted.WeightedSum.shortLabel,...    
%     combiners.DlrCombiner.shortLabel,...
%     combiners.bcc.IbccSampling.shortLabel,...
%     combiners.SimpleMajorityVoting.shortLabel, ...
%     combiners.weighted.WeightedMajority.subShortLabel, ...
%     combiners.bcc.IbccVbEmInit.shortLabel, ...
%     combiners.bcc.IbccEm.shortLabel, ...    
    };

runner = ExpRunner(expSettings);

%arg1 = new data for base classifiers?
%arg2 = whether to run agents again, load saved, run fake agents, do
%nothing
%arg3 = draw graphs?
%arg4 = list of combination methods to use
%arg5 = sort results: use mean to sort results so that those with similar
%value are placed next to one another in the graph. Set to false to
%preserve the real order of data points as presented to the combiners.
runner.Alpha0 = [0.2 0.2; 0.2 0.2];
[combinedPost, agentRatings, labels] = runner.runCombineAll(true, 'newFakeAgents', true, combMethods, false);
% [combinedPost, agentRatings, labels] = runner.runCombineAll(false, 'loadSaved', true, combMethods, false);
graphs.ClassifierPerformanceGraph.drawRoc(combinedPost, labels(1:expSettings.nSamples)-1, combMethods, false);