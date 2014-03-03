%Run this code from this directory
%clear
import settings.*

%load('/homes/49/edwin/results/agent_clustering_meam-magnitude_ordered_relabelling/data/5.mat');

%experimentalSettingsCombination
%bccSettingsInversions;
% bccSettings;
bccSettingsDynamicsExp3;
expSettings.topSaveDir = [expSettings.topSaveDir '/singletest'];
expSettings.definiteSwitch = [0.07*5 0.08*5 0.09*5 0.10*5 0.13*5 0.14*5 0.15*5 0.16*5];
expSettings.iPropKnown = 3;

%cluster agents over a data sequence of 20. use posteriors as cluster input
%data.
% expSettings.clusterData = 'post';
% expSettings.weightClusterInput = false;
% expSettings.topSaveDir = '/homes/49/edwin/results/combination/inf5_noninf11_agents5dlr';
% expSettings.expLabel = 'dlrCombiners';
%reuseSubs = true;

%pick the combination methods to use.
combMethods = {...
     combiners.MeanDecision.shortLabel,...    
     combiners.bcc.DynIbccVb.shortLabel,...
     combiners.bcc.IbccVb.shortLabel,...
...%      combiners.bcc.IbccEm.shortLabel,...
...%      combiners.SimpleMajorityVoting.shortLabel, ...     
...%      combiners.weighted.WeightedSum.shortLabel,...
     combiners.DlrCombiner.shortLabel ...
...%      combiners.weighted.NaiveBayesLogodds.precalcShortLabel, ...
%      combiners.weighted.WeightedSum.shortLabel,...
%      combiners.weighted.WeightedMajority.subShortLabel...
%      combiners.bcc.Ibcc.shortLabel, ...
    };

runner = ExpRunner(expSettings);
[combinedPost, agentRatings, labels] = runner.runCombineAll(true, 'newAgents', true, combMethods);
% [combinedPost, agentRatings, labels] = runner.runCombineAll(false, 'loadSaved', true, combMethods);

%set newdata argument to false to reuse the previous data; newagents
%argument to false to reuse existing agent-sensor assignments and training.
%runner.runCombineByCluster(true, true, combMethods);
%runner.runCombineByCluster(false, false, combMethods);
graphs.ClassifierPerformanceGraph.drawRoc(combinedPost, labels(1:expSettings.nSamples)-1, combMethods, false);