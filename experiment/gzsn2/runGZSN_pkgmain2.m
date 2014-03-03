function runGZSN_pkgmain2(rootDir)

display('starting GZSN classifier combination');

%If the paths are added to the default MATLAB path we can comment these out
%and save a bit of time
if nargin<1
    rootDir = '/homes/49/edwin/matlab/distrib/gzsnBundle/';
    addpath([rootDir '../../../matlab']);
end

% !!! Any settings from the old package main that we haven't passed through
% should be put in a new settings file here
settings.gzsn_pkg;

%pick the combination methods to use.
combMethods = {...
...%    combiners.MeanDecision.shortLabel,...  
...%     combiners.SimpleMajorityVoting.shortLabel, ...
...%     combiners.weighted.WeightedMajority.subShortLabel, ...
     combiners.bcc.IbccVb.shortLabel,...   
    };

interpretEmptyScreenLabels = false;
legalScores = [-1 1 3];
runExperiment

end
