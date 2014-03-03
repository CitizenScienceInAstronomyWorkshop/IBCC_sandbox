%TRY ro reach 95% accuracy without adding in the assumption about
%unlabelled objects. Try to get good performance by doing subsampling that
%would be possible without knowing the test labels a priori, e.g. using small
%number of test samples against lots of training data and balancing the
%training data. So just the test portion is random. Or e.g. filtering
%agents with few responses.

if ~exist('rootDir','var')
    rootDir = [];%'/home/edwin/matlab/data/galaxyZoo3/';
end
expSet = settings.ExpSettings(rootDir, 'gzsnData.csv');

expSet.nRepeats = 1;
expSet.topSaveDir = '/home/edwin/matlab/data/galaxyZoo3/combiner_output/';%exp2its';%
%expSet.topSaveDir = '/homes/49/edwin/matlab/data/revisedVB/exp3';
%expSet.noOverwriteData = true;
%expSet.expLabel = 'dlrCombiners'; %experiment label
expSet.expLabel = 'bcc';
expSet.nScores = 3;
expSet.nFolds = 10;

expSet.outputDir = '/home/edwin/matlab/data/galaxyZoo3/combiner_output/';
expSet.saveAs = '';

expSet.includeUnlabelled = true;
expSet.batches = false;

%if we exclude unlabelled but do no further subsampling, IBCCVB reaches 76%
%AUC. 

%filter out data with not enough responses from frequently-responding agents
expSet.minAgentResp = 0; %min no. responses for a frequently-responding agent 
expSet.minFreqAgents = 0; %min number of freqently-responding agents per data point
expSet.startNegAssets = 1;
expSet.maxNoAssets = 500;

bccSet = settings.BccSettings();
bccSet.screenLabelsAsAgent = false; %*****

% bccSet.Alpha = [2 1.5 1; 1 1.5 2];

% bccSet.Alpha = [20 1 1; 10 1 2];
% bccSet.Alpha = [0.6 0.5 0.5; 0.5 0.5 0.6];
% bccSet.Alpha = [0.55 0.4 0.05; 0.2 0.4 0.4];
% bccSet.Alpha = [0.5 0.36 0.09; 0.24 0.36 0.35];
bccSet.Alpha = [0.5 0.38 0.12; 0.36 0.3 0.34];
% bccSet.Alpha = [0.5 0.3 0.05; 0.18 0.36 0.41]; %best so far on no subsampling(auc=0.9042 with IbccDiff) %standard priors: %[0.66 0.5 0.4; 0.49 0.43 0.44].*50; [0.6 0.5 0.5; 0.48 0.49 0.59].*50;
bccSet.trustedAlpha = [0.6 0.3 0.1; 0.2 0.3 0.5].*10;%[0.6 0.5 0.45; 0.45 0.5 0.6].*10;%0.48 0.49 0.59].*10;

bccSet.Alpha(:,:) = bccSet.Alpha(:,:) .* 10;
bccSet.Alpha(:,1) = bccSet.Alpha(:,1) + 50;

%use this for the "real" class proportions
bccSet.nu = {[100 100]}; %***? keep this??

bccSet.minScore = -1;
bccSet.maxScore = 3;
bccSet.maxIt = 500; %***
bccSet.convThreshold = 10^-4;
bccSet.priorAdjLevel = 0;%75;

%IDEALLY: turn off pooling in IBCC diff. Subsample without screening labels
%gives good normal IBCCVB results. No subsampling without screening labels
%gives >80% results with IBCC diff. Otherwise, use screening labels, but
%mean needs to be >10% worse.

% ***the line below, and +1s to the alphas
% bccSet.AlphaDiff = ([0.6 0.5 0.1; 0.6 0.5 0.1; 0.48 0.49 0.59; 0.48 0.49 0.59].*1); %decreased p(least likely response)
% bccSet.AlphaDiff = [0.5 0.3 0.05;0.5 0.3 0.05; 0.18 0.36 0.41;0.18 0.36 0.41];
bccSet.AlphaDiff = ([0.6 0.5 0.5; 0.6 0.5 0.5; 0.5 0.5 0.6; 0.5 0.5 0.6].*1);%previously 10 
bccSet.AlphaDiff(2,:)=bccSet.AlphaDiff(2,:).*0.1; 
bccSet.AlphaDiff(4,:)=bccSet.AlphaDiff(4,:).*0.1;

bccSet.IbccMapFunction = @mapGZSNScoresToVotes;

bccSet.debug = false;

filterPos = true;
interpretEmptyScreenLabels = false;
keepInfreqAgents = true;