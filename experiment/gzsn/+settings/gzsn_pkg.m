if ~exist('rootDir','var')
    rootDir = [];
end
expSet = settings.ExpSettings(rootDir, 'input.csv', 'data/gzsnCombinedScores');

expSet.drawGraphs = false;
expSet.writeTrainingResults = true;

expSet.nRepeats = 1;
expSet.expLabel = 'bcc';
expSet.nScores = 3;
expSet.nFolds = 1;

expSet.saveAs = 'csv';

expSet.includeUnlabelled = true;
expSet.batches = false;

%filter out data with not enough responses from frequently-responding agents
expSet.minAgentResp = 1; %min no. responses for a frequently-responding agent 
expSet.minFreqAgents = 1; %min number of freqently-responding agents per data point
expSet.startNegAssets = 1;
expSet.maxNoAssets = 0;

bccSet = settings.BccSettings();
bccSet.screenLabelsAsAgent = true;
bccSet.Alpha = [0.5 0.3 0.05; 0.18 0.36 0.41]; %best so far on no subsampling with original dataset including interpret empty screening labels as negative examples(auc=0.9042 with IbccDiff) %standard priors: %[0.66 0.5 0.4; 0.49 0.43 0.44].*50; [0.6 0.5 0.5; 0.48 0.49 0.59].*50;
bccSet.trustedAlpha = [0.7 0.5 0.5; 0.48 0.49 0.69].*20;

%use this for the "real" class proportions
% expSet.nu = {[40000 4]};
bccSet.nu = {[100 100]};

bccSet.minScore = -1;
bccSet.maxScore = 3;
bccSet.maxIt = 500;

bccSet.AlphaDiff = ([0.6 0.5 0.5; 0.6 0.5 0.5; 0.48 0.49 0.59; 0.48 0.49 0.59].*50); 
% bccSet.AlphaDiff(2,:)=bccSet.AlphaDiff(2,:).*0.7; 
bccSet.AlphaDiff(4,:)=bccSet.AlphaDiff(4,:).*0.025;