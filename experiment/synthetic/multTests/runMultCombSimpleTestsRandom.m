%Run this code from this directory
%clear
import settings.*

%load('/homes/49/edwin/results/agent_clustering_meam-magnitude_ordered_relabelling/data/5.mat');

bccSettingsSimpleRandom;

if ~exist(expSettings.getCombinerDir(), 'file')
    mkdir(expSettings.getCombinerDir());
end

%pick the combination methods to use.
combMethods = {...
     combiners.MeanDecision.shortLabel,...
     combiners.bcc.Ibcc.shortLabel,...      
%           combiners.bcc.IbccVb.shortLabel,...     
%      combiners.bcc.IbccVbAux.shortLabel,...     
%      combiners.weighted.WeightedMajority.subShortLabel...
%      combiners.weighted.WeightedSum.shortLabel,...     
%      combiners.bcc.IbccVbEmInit.shortLabel, ...
%      combiners.bcc.IbccEm.shortLabel, ...     
%     combiners.SimpleMajorityVoting.shortLabel, ...
    };

Alpha0 = [10 1; 1 10];

multCombSimpleTests;