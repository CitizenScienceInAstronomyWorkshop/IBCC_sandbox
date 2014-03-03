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

%filter out data with not enough responses from frequently-responding agents
expSet.minAgentResp = 0; %min no. responses for a frequently-responding agent 
expSet.minFreqAgents = 0; %min number of freqently-responding agents per data point
expSet.startNegAssets = 1;
expSet.maxNoAssets = 0;

expSet.recreateFolds = true;

bccSet = settings.BccSettings();
bccSet.screenLabelsAsAgent = false; %*****

% bccSet.Alpha = [2 1.5 1; 1 1.5 2];

bccSet.Alpha = [0.6 0.5 0.3; 0.35 0.45 0.6]+ 1;% .* 0.5;
% bccSet.Alpha = [18.5 0.3 0.2; 18.3 0.3 0.3] .* 10;

% bccSet.Alpha = [0.55 0.4 0.05; 0.2 0.4 0.4];
% bccSet.Alpha = [0.5 0.3 0.15; 0.18 0.36 0.40]; %best so far on no subsampling(auc=0.9042 with IbccDiff) %standard priors: %[0.66 0.5 0.4; 0.49 0.43 0.44].*50; [0.6 0.5 0.5; 0.48 0.49 0.59].*50;
bccSet.trustedAlpha = [0.6 0.5 0.45; 0.45 0.5 0.6];%0.48 0.49 0.59].*10;

% bccSet.Alpha = [1 1 1; 1 1 1];
% bccSet.Alpha = [0.5 0.3 0.2; 0.3 0.3 0.3];
% bccSet.Alpha = bccSet.Alpha  .* 10;
% bccSet.Alpha(:,1) = bccSet.Alpha(:,1) + 60;%1;

%use this for the "real" class proportions
bccSet.nu = {[100 500]}; %***? keep this??

bccSet.minScore = -1;
bccSet.maxScore = 3;
bccSet.maxIt = 200; %***
bccSet.fixedNIt = 500;
bccSet.convThreshold = 10^-6;
% bccSet.flipSamples = false;

%IDEALLY: turn off pooling in IBCC diff. Subsample without screening labels
%gives good normal IBCCVB results. No subsampling without screening labels
%gives >80% results with IBCC diff. Otherwise, use screening labels, but
%mean needs to be >10% worse.

% ***the line below, and +1s to the alphas
% bccSet.AlphaDiff = ([0.6 0.5 0.1; 0.6 0.5 0.1; 0.48 0.49 0.59; 0.48 0.49 0.59].*1); %decreased p(least likely response)
% bccSet.AlphaDiff = [0.5 0.3 0.05;0.5 0.3 0.05; 0.18 0.36 0.41;0.18 0.36 0.41];
% bccSet.AlphaDiff = ([1.6 1.5 1.3; 0.6 0.5 0.3; 1.3 1.5 1.6; 0.3 0.5 0.6].*1);%previously 10 
bccSet.AlphaDiff = ([1.6 1.5 1.3; 1.4 1.6 1.7; 1.3 1.6 1.8; 1 1 1].*1);%

% bccSet.AlphaDiff(1,:)=bccSet.AlphaDiff(1,:).*10 + 20; 
% bccSet.AlphaDiff(3,:)=bccSet.AlphaDiff(3,:).*10 + 20;

% bccSet.AlphaDiff(2,:)=bccSet.AlphaDiff(2,:).*0.1; 
bccSet.AlphaDiff(4,:)=bccSet.AlphaDiff(4,:).*0.1;

bccSet.IbccMapFunction = @mapGZSNScoresToVotes;

bccSet.debug = false;