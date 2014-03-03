
%bcc multiple test settings: easy tests

expSettings = settings.ExpSettings();

%Ibcc
expSettings.propKnown = [0.1 0.2 0.5 0.8 1];%[0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45];%[0 0.05 0.15 0.25];
%expSettings.propKnown = [0 0.1 0.4];
%expSettings.propKnown = [0 0.1 0.25 0.4];

expSettings.iPropKnown = 2;

expSettings.lambdaSym = [0.1, 0.2, 0.3, 0.4, 0.49, 0.5, 0.6, 0.7, 0.8, 0.9];
expSettings.iLambdaSym = 4;

expSettings.lambdaMag = [1 10 50 100];
expSettings.iLambdaMag = 2;

expSettings.nu = {[10 10]};
expSettings.iNu = 1;

%detailed runs use these settings
%propKnown = [0 0.005 0.01 0.02 0.04 0.06 0.08 0.1 0.2];
%lambdaSym = [0.5, 0.6, 0.7, 0.8, 0.9];
%lambdaMag = [0.1 0.5 1 5];

expSettings.nDatasets = 5;
expSettings.nRepeats = 1;

expSettings.nSamples = 500;

expSettings.initCombinedPosteriors();

%cluster agents over a data sequence of 20. use posteriors as cluster input
%data.
expSettings.clusterData = 'post';
expSettings.weightClusterInput = false;
expSettings.topSaveDir = '/homes/49/edwin/matlab/combination/data/dynamicVB/simple/exp3dyn_playground';
%expSettings.noOverwriteData = true;
%expSettings.expLabel = 'dlrCombiners'; %experiment label
expSettings.expLabel = 'combiners';

%other more boring things
expSettings.saveImgs = true;
expSettings.agentType='dlr';

%data set sizes
expSettings.nTrainingSamples = 0;
expSettings.lengthHistory = 50;

%data generation
expSettings.means = [-3 -2 -1 ; 3 2 1 ];

%test with sensors all of same quality
expSettings.deviations = [1 2 1 ; 1 1 2];

expSettings.p_c1 = 0.5;

expSettings.noninfMean = 0;
expSettings.noninfDev = 5;

%sensors

%settings used: 
% inf =5 noninf = 11; max =10 min = 3
% inf = 3 noninf = 8
% inf = 2, noninf = 2; max = 3, min = 1
% inf = 1 noninf = 1

expSettings.nInfSensors = 3;
expSettings.nNoninfSensors = 27;

expSettings.pSwitch = 0.0; % probability that a sensor will stop being useful
expSettings.definiteSwitch = [0.33 0.67];%[0.19 0.39 0.59 0.79];
display('Changepoints are not saved, so the data should be re-generated if you want to plot them');
expSettings.pMissing = 0.00;
expSettings.pCorrupt = 0.0;
expSettings.pFlipLabel = 0.0;
expSettings.pFlipSensor = 0.0;
expSettings.spreadTargets = true;

%agents and sensor allocation
expSettings.nAgents = 40;
expSettings.nInformedAgents = 4;

expSettings.maxSensors = 6;
expSettings.minSensors = 2;

%number of clusters and clustering window length
expSettings.seqLength = 20;
expSettings.K = 2;

expSettings.fakeAgents = zeros(2, 2, expSettings.nAgents);
expSettings.fakeAgents(:, :, 1) = [ 0.9 0.1; 0.1 0.9]; % 90% correct agent
expSettings.fakeAgents(:, :, 2) = [ 0.1 0.9; 0.1 0.9 ]; % 90% correct agent
expSettings.fakeAgents(:, :, 3) = [ 0.01 0.99; 0.19, 0.81 ]; % 90% correct agent
expSettings.fakeAgents(:, :, 4) = [ 0.01 0.99; 0.19, 0.81 ]; % 90% correct agent

for a=expSettings.nInformedAgents+1:expSettings.nAgents
    expSettings.fakeAgents(:, :, a) = [ 0.5 0.5; 0.5 0.5 ]; % useless agent
end
% 
% import combiners.ibccsampling.GibbsSampling
% GibbsSampling.gibbsSamples = 100;
% GibbsSampling.sampleInterval = 1;
% GibbsSampling.burnIn = 30;