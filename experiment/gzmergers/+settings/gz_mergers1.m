rootDir = '/homes/49/edwin/matlab/data/galaxyZooMergers/';

expSet = settings.ExpSettings(rootDir);

expSet.nRepeats = 1;
expSet.topSaveDir = [rootDir 'combiner'];%exp2its';%

%expSet.dontOverwrite = true;
expSet.expLabel = 'ibcc';
expSet.nFolds = 5;
expSet.nScores = 3;
expSet.includeUnlabelled = false;
expSet.batches = false;
expSet.writeTrainingResults = true;

expSet.saveAs = '';

%filter out data with not enough responses from frequently-responding agents
expSet.minAgentResp = 1; %min no. responses for a frequently-responding agent 
expSet.minFreqAgents = 1; %min number of freqently-responding agents per data point
expSet.startNegAssets = 1;
expSet.maxNoAssets = 0;

bccSet = settings.BccSettings();
bccSet.convThreshold = 10^-3;
%use this for the "real" class proportions
% expSet.nu = {[40000 4]};
bccSet.nu = {[3000 600]};

bccSet.IbccMapFunction = @mapGZMergerScores; %for using all 3 scores
bccSet.scoreMap = [1 2 3; 1 2 3];
bccSet.minScore = 1;
bccSet.maxScore = 3;
bccSet.maxIt = 50;

bccSet.screenLabelsAsAgent = false;
% bccSet.Alpha = [1.1 1 0.9; 0.9 1 1.1] .* 1000; 
% bccSet.Alpha = [0.6 0.5 0.5; 0.48 0.49 0.59].*50;
% bccSet.Alpha = [7 6 2; 4 5 6].*10;
bccSet.Alpha = ([19 3.5 0.5; 17.5 3.0 2.5] .* 10); %([19 3.5 0.5; 17.5 3.0 2.5] .* 10);
% bccSet.AlphaDiff = ([19.4 3.2 0.4; 16.5 3.4 3.4; 16.5 3.4 3.4; 18.2 3.3 1.5] .* 5); %([19 3.5 0.5; 17.5 3.0 2.5] .* 10);
bccSet.AlphaDiff = ([19 3.5 0.5; 19 3.5 0.5; 17.5 3.0 2.5; 17.5 3.0 2.5].*20); 
% bccSet.AlphaDiff(2,:)=bccSet.AlphaDiff(2,:).*0.025; 
bccSet.AlphaDiff(4,:)=bccSet.AlphaDiff(4,:).*0.025;
bccSet.priorAdjLevel = 50;%100; results in the email use this setting. We ran a test changing this value.

% bccSet.Alpha = [1 1 1; 1 1 1] .* 50;
bccSet.trustedAlpha = [0.6 0.5 0.5; 0.48 0.49 0.59].*50;
bccSet.useLikelihood = false;
