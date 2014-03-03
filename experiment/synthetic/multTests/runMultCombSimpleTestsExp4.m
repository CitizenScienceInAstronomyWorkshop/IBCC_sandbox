%Run this code from this directory
%clear
import settings.*

%load('/homes/49/edwin/results/agent_clustering_meam-magnitude_ordered_relabelling/data/5.mat');

bccSettingsSimple4;

if ~exist(expSettings.getCombinerDir(), 'file')
    mkdir(expSettings.getCombinerDir());
end

%pick the combination methods to use.
combMethods = {...
     combiners.MeanDecision.shortLabel,...
     combiners.DlrCombiner.shortLabel, ...
     combiners.bcc.IbccVb.shortLabel,...
     combiners.bcc.IbccEm.shortLabel, ...
     combiners.bcc.Ibcc.shortLabel, ...
     combiners.SimpleMajorityVoting.shortLabel, ...
     combiners.weighted.WeightedSum.shortLabel,...
     combiners.weighted.WeightedMajority.subShortLabel...
    };

multCombSimpleTests;