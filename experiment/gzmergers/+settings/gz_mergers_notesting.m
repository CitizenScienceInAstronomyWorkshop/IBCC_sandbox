expSet = settings.ExpSettings('/homes/49/edwin/matlab/data/galaxyZooMergers/', [], outputDirLabel);

expSet.nRepeats = 1;
expSet.topSaveDir = '/homes/49/edwin/matlab/data/galaxyZooMergers/combiner';%exp2its';%

%expSet.dontOverwrite = true;
expSet.expLabel = 'ibcc';
expSet.nFolds = 1;
expSet.nScores = 3;
expSet.includeUnlabelled = true;
expSet.batches = false;
expSet.writeTrainingResults = true;

expSet.saveAs = 'csv';

%filter out data with not enough responses from frequently-responding agents
expSet.minAgentResp = 1; %min no. responses for a frequently-responding agent 
expSet.minFreqAgents = 1; %min number of freqently-responding agents per data point
expSet.startNegAssets = 1;
expSet.maxNoAssets = 0;

bccSet = settings.BccSettings();
bccSet.convThreshold = 10^-4;
%use this for the "real" class proportions
% expSet.nu = {[40000 4]};
bccSet.nu = {[1000 1000]};

bccSet.IbccMapFunction = @mapGZMergerScores; %for using all 3 scores
bccSet.scoreMap = [1 2 3; 1 2 3];
bccSet.minScore = 1;
bccSet.maxScore = 3;
bccSet.maxIt = 500;

bccSet.screenLabelsAsAgent = false;
% bccSet.Alpha = [1.1 1 0.9; 0.9 1 1.1] .* 1000; 
% bccSet.Alpha = [0.6 0.5 0.5; 0.48 0.49 0.59].*50;
% bccSet.Alpha = [7 6 2; 4 5 6].*10;
bccSet.Alpha = ([19 3.5 0.5; 17.5 3.0 2.5] .* 0.3);
% bccSet.Alpha = [1 1 1; 1 1 1] .* 50;
bccSet.trustedAlpha = [0.6 0.5 0.5; 0.48 0.49 0.59].*50;
bccSet.useLikelihood = false;
