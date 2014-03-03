expSettings = settings.ExpSettings();

expSettings.flipAgents = 0.5;

expSettings.nDatasets =25;
expSettings.nRepeats = 1;

expSettings.nSamples = 50;

%Ibcc
expSettings.propKnown = [0 0.1 0.2];%[0 0.05 0.1 0.15 0.2 0.25];
expSettings.iPropKnown = 1;

%the lower this is, the more we reinforce the belief that the classifiers
%are 'good'.
expSettings.lambdaSym = [0.1, 0.2, 0.3, 0.4, 0.49, 0.5, 0.6, 0.7, 0.8, 0.9];
expSettings.iLambdaSym = 5;

%the lower this is, the more weight we give to the priors. For VB, this
%also gives more weight to the data, i.e. should converge quicker.
expSettings.lambdaMag = [0.001 0.1 0.5 1 5 10 50 100];
expSettings.lambdaMag = expSettings.lambdaMag .* 100 ./ expSettings.nSamples;
expSettings.iLambdaMag = 4;

%expSettings.nu = {[20 20]};

expSettings.initCombinedPosteriors();

expSettings.topSaveDir = '/homes/49/edwin/matlab/combination/data/vbIbccPaper/mini';%exp2_inv';%exp2inv_itstest';%exp2_inv_nknownsmall';
%expSettings.topSaveDir = '/homes/49/edwin/matlab/combination/data/revisedVB/exp3';
expSettings.noOverwriteData = true;
%expSettings.expLabel = 'dlrCombiners'; %experiment label
expSettings.expLabel = 'bcc';

%cluster agents over a data sequence of 20. use posteriors as cluster input
%data.
expSettings.clusterData = 'post';
expSettings.weightClusterInput = false;

%other more boring things
expSettings.saveImgs = true;
expSettings.agentType='dlr';

expSettings.nTrainingSamples = 20;
expSettings.lengthHistory = 50;

%data generation
expSettings.means = [-1 -1 -1 -1 -1 -1 -1 -1 -1 -1; 1 1 1 1 1 1 1 1 1 1];

%test with sensors all of same quality
expSettings.deviations = [1 1 1 1 1 1 1 1 1 1; 2 2 2 2 2 2 2 2 2 2];

expSettings.p_c1 = 0.5;

expSettings.noninfMean = 0;
expSettings.noninfDev = 10;

%sensors

%settings used: 
% inf =5 noninf = 11; max =10 min = 3
% inf = 3 noninf = 8
% inf = 2, noninf = 2; max = 3, min = 1
% inf = 1 noninf = 1

expSettings.nInfSensors = 20; %10
expSettings.nNoninfSensors = 40; %20

expSettings.pSwitch = 0.0; % probability that a sensor will stop being useful
expSettings.pMissing = 0.00;
expSettings.pCorrupt = 0.0;
expSettings.pFlipLabel = 0.0;
expSettings.pFlipSensor = 0.0;

%agents and sensor allocation
expSettings.nAgents = 40;
expSettings.nInformedAgents = 16;

expSettings.maxSensors = 6;
expSettings.minSensors = 2;

%number of clusters and clustering window length
expSettings.seqLength = 20;
expSettings.K = 2;