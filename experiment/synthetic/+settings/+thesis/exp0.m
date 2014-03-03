
%This is the basic template for the experiments. Start from this, and then
%alter whichever independent variable belongs to the current experiment.

% expSet.initCombinedPosteriors();

rootForAllExpts = '/homes/49/edwin/matlab/data/thesis_synth_sharing/';

expSet = settings.ExpSettings([rootForAllExpts 'exp0/']);

expSet.noScore = -1;

expSet.propKnown = 0.025;
expSet.iPropKnown = 1;

expSet.nRepeats = 1;

expSet.nSamples = 1000; %******

expSet.topSaveDir = expSet.rootDir;

expSet.expLabel = 'th1';

%other more boring things
expSet.saveImgs = true;

%Settings for dataset synthesisation
synthSet = settings.SynthSettings(expSet.topSaveDir);

%--- No. DATASETS --------------------------------------------------------------
synthSet.nDatasets = 25;  %******
%-------------------------------------------------------------------------------
display('set to overwrite data');
synthSet.noOverwriteData = false;

synthSet.agentType='logistic';

%data set sizes
synthSet.nTrainingSamples = 5;

display('Think we need to change this: want to have fixed agents and change the sensors instead');
synthSet.lengthHistory = 100;

synthSet.nAgents = 5;

%data generation
synthSet.means = ones(2, synthSet.nAgents);
synthSet.means(2,:) = -synthSet.means(1,:);%[1 1 1 1 1; -1 -1 -1 -1 -1];

%test with sensors all of same quality
synthSet.deviations = 1 .* ones(2, synthSet.nAgents);
synthSet.deviations(2,:) = -synthSet.deviations(1,:);%[1 1.5 3 9 10; 1 1.5 3 9 10];

synthSet.cal = ones(1,synthSet.nAgents);%[5 5 2.5 0.5 0.5];
synthSet.bias = zeros(1, synthSet.nAgents);%[0 0 0 0 0];

synthSet.p_c1 = 0.5;

synthSet.noninfMean = 0;
synthSet.noninfDev = 5;

%sensors
synthSet.nInfSensors = 5;
synthSet.nNoninfSensors = 0;

synthSet.pSwitch = 0.0; % probability that a sensor will stop being useful
synthSet.pMissing = 0.00;
synthSet.pCorrupt = 0.0;
synthSet.pFlipLabel = 0.0;
synthSet.pFlipSensor = 0.0;

%agents and sensor allocation
synthSet.nInformedAgents = synthSet.nAgents;

synthSet.maxSensors = 1;
synthSet.minSensors = 1;

bccSet = settings.BccSettings();

bccSet.nu = {[100 100]};
bccSet.iNu = 1;

bccSet.Alpha = [12 8; 8 12];