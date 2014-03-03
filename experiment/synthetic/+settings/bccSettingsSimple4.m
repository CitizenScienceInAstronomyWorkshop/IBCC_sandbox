
%bcc multiple test settings: easy tests

expSettings = settings.ExpSettings();

%Ibcc
expSettings.propKnown = [0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45];%[0 0.05 0.15 0.25];
%expSettings.propKnown = [0 0.1 0.4];
%expSettings.propKnown = [0 0.1 0.25 0.4];

expSettings.iPropKnown = 2;

expSettings.lambdaSym = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
expSettings.iLambdaSym = 2;

expSettings.lambdaMag = [1 10 50 100];
expSettings.iLambdaMag = 2;

%detailed runs use these settings
%propKnown = [0 0.005 0.01 0.02 0.04 0.06 0.08 0.1 0.2];
%lambdaSym = [0.5, 0.6, 0.7, 0.8, 0.9];
%lambdaMag = [0.1 0.5 1 5];

expSettings.nDatasets = 10;
expSettings.nRepeats = 3;

expSettings.nSamples = 500;

expSettings.initCombinedPosteriors();

%cluster agents over a data sequence of 20. use posteriors as cluster input
%data.
expSettings.clusterData = 'post';
expSettings.weightClusterInput = false;
expSettings.topSaveDir = '/homes/49/edwin/matlab/combination/data/simple_base_classifiers/exp4';
expSettings.noOverwriteData = true;
%expSettings.expLabel = 'dlrCombiners'; %experiment label
expSettings.expLabel = 'combiners';

%other more boring things
expSettings.saveImgs = true;
expSettings.agentType='dlr';

%data set sizes
expSettings.nTrainingSamples = 0;
expSettings.lengthHistory = 50;

%data generation
expSettings.means = [-3 -2 -1 3 2 -1 -1 -0.5 -0.5 -1 -3 -2 -1 3 2 -1 -1 -0.5 -0.5 -1;...
    3 2 1 -3 -2 1 1 0.5 0.5 1 3 2 1 -3 -2 1 1 0.5 0.5 1];

%test with sensors all of same quality
expSettings.deviations = [1 2 1 3 1 2 2 2 2 2 1 2 1 3 1 2 2 2 2 2; ...
    1 1 2 1 3 2 2 2 2 2 1 1 2 1 3 2 2 2 2 2];

expSettings.p_c1 = 0.5;

expSettings.noninfMean = 0;
expSettings.noninfDev = 5;

%sensors
expSettings.nInfSensors = 20;
expSettings.nNoninfSensors = 0;

expSettings.pSwitch = 0.0; % probability that a sensor will stop being useful
expSettings.pMissing = 0.00;
expSettings.pCorrupt = 0.0;
expSettings.pFlipLabel = 0.0;
expSettings.pFlipSensor = 0.0;

%agents and sensor allocation
expSettings.nAgents = 40;
expSettings.nInformedAgents = 4;

expSettings.maxSensors = 6;
expSettings.minSensors = 2;

%number of clusters and clustering window length
expSettings.seqLength = 20;
expSettings.K = 2;

expSettings.fakeAgents = zeros(2, 2, expSettings.nAgents);
expSettings.fakeAgents(:, :, 1) = [ 0.9 0.1; 0.1 0.9 ]; % 90% correct agent
expSettings.fakeAgents(:, :, 2) = [ 0.9 0.1; 0.1 0.9 ]; % 90% correct agent
expSettings.fakeAgents(:, :, 3) = [ 0.99 0.01; 0.19, 0.81 ]; % 90% correct agent
expSettings.fakeAgents(:, :, 4) = [ 0.8 0.2; 0.2, 0.8 ]; % 90% correct agent
expSettings.fakeAgents(:, :, 5) = [ 0.7 0.3; 0.3, 0.7 ]; % 90% correct agent
expSettings.fakeAgents(:, :, 6) = [ 0.6 0.4; 0.4, 0.6 ]; % 90% correct agent
expSettings.fakeAgents(:, :, 7) = [ 0.5 0.5; 0.2, 0.8 ]; % 90% correct agent
expSettings.fakeAgents(:, :, 8) = [ 0.8 0.2; 0.5, 0.5 ]; % 90% correct agent
expSettings.fakeAgents(:, :, 9) = [ 0.75 0.25; 0.75, 0.25 ]; % 90% correct agent
expSettings.fakeAgents(:, :, 10) = [ 0.85 0.15; 0.85, 0.15 ]; % 90% correct agent

for a=11:20
    expSettings.fakeAgents(:, :, a) = flipdim(expSettings.fakeAgents(:, :, a-10), 2); % useless agent
end

import combiners.ibccsampling.GibbsSampling
GibbsSampling.gibbsSamples = 100;
GibbsSampling.sampleInterval = 1;
GibbsSampling.burnIn = 30;