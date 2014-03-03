expSet = settings.ExpSettings('/homes/49/edwin/matlab/data/zach/', '/datacsv.csv');

expSet.nRepeats = 1;
expSet.topSaveDir = '/homes/49/edwin/matlab/data/zach/combination';%exp2its';%

%expSet.dontOverwrite = true;
expSet.expLabel = 'ibcc';
expSet.nFolds = 5;
expSet.nScores = 2;
expSet.includeUnlabelled = true;
expSet.batches = false;
expSet.writeTrainingResults = true;

expSet.saveAs = 'csv';

%filter out data with not enough responses from frequently-responding agents
expSet.minAgentResp = 0; %min no. responses for a frequently-responding agent 
expSet.minFreqAgents = 0; %min number of freqently-responding agents per data point
expSet.startNegAssets = 1;
expSet.maxNoAssets = 0;

bccSet = settings.BccSettings();
bccSet.screenLabelsAsAgent = false;

bccSet.IbccMapFunction = @mapZachScores; %for using all 3 scores

% bccSet.scoreMap = [1 2; -1 1];
bccSet.minScore = -1;
bccSet.maxScore = 1;

bccSet.maxIt = 100;
bccSet.convThreshold = 10^-3;

bccSet.nu = {[1000 1000]};
bccSet.Alpha = ([4 12; 4 12] .* 0.1);
bccSet.priorAdjLevel = 0;