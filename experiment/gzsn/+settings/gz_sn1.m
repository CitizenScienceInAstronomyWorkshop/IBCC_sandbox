if ~exist('rootDir','var')
    rootDir = [];
end
expSet = settings.ExpSettings(rootDir, 'gzsnData.csv');

expSet.nRepeats = 1;
expSet.topSaveDir = '/homes/49/edwin/matlab/data/galaxyZoo3/combiner_output';%exp2its';%
%expSet.topSaveDir = '/homes/49/edwin/matlab/data/revisedVB/exp3';
%expSet.noOverwriteData = true;
%expSet.expLabel = 'dlrCombiners'; %experiment label
expSet.expLabel = 'bcc';
expSet.nScores = 3;
expSet.nFolds = 10;

expSet.outputDir = '/homes/49/edwin/matlab/data/galaxyZoo3/combiner_output';
expSet.saveAs = '';

expSet.includeUnlabelled = true;
expSet.batches = false;

%filter out data with not enough responses from frequently-responding agents
expSet.minAgentResp = 0; %min no. responses for a frequently-responding agent 
expSet.minFreqAgents = 0; %min number of freqently-responding agents per data point
expSet.startNegAssets = 1;
expSet.maxNoAssets = 0;

bccSet = settings.BccSettings();
bccSet.screenLabelsAsAgent = false;%true

% bccSet.Alpha = [2 1.5 1; 1 1.5 2];
% bccSet.Alpha = ([0.6 0.5 0.5; 0.48 0.49 0.59].*1); 

bccSet.Alpha = [0.5 0.3 0.05; 0.18 0.36 0.41];% + 1; %added 1 to current best %best so far on no subsampling(auc=0.9042 with IbccDiff) %standard priors: %[0.66 0.5 0.4; 0.49 0.43 0.44].*50; [0.6 0.5 0.5; 0.48 0.49 0.59].*50;
bccSet.trustedAlpha = [0.6 0.5 0.5; 0.48 0.49 0.59].*10; %(increased trust in screening step)

%use this for the "real" class proportions
bccSet.nu = {[100 10]};

bccSet.minScore = -1;
bccSet.maxScore = 3;
bccSet.maxIt = 100;

% bccSet.AlphaDiff = ([0.6 0.5 0.1; 0.6 0.5 0.1; 0.48 0.49 0.59; 0.48 0.49 0.59].*1); %decreased p(least likely response)
bccSet.AlphaDiff = ([0.6 0.5 0.5; 0.6 0.5 0.5; 0.5 0.5 0.6; 0.5 0.5 0.6].*1); 
bccSet.AlphaDiff(2,:)=bccSet.AlphaDiff(2,:).*0.7; 
bccSet.AlphaDiff(4,:)=bccSet.AlphaDiff(4,:).*0.025;

bccSet.IbccMapFunction = @mapGZSNScoresToVotes;

bccSet.debug = false;