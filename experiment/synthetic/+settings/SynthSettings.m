classdef SynthSettings < handle
    %SYNTHSETTINGS Settings for experiments using synthetic data: settings
    %to generate the data
    
    properties
                
        nDatasets = 10;        
        
        %DATA DIR/FILENAMES
        topSaveDir = 'unknownSaveDirectory';
        multCombTestFile = 'multCombinerTests.mat';
        
        %AGENTS
        agentType = 'dlr';
                
        nAgents = 5;
        nInformedAgents = -1;
        
        fakeAgents = [];
        
        %proportion of agents whose outputs will be flipped
        flipAgents = 0;
        
        %size of data window for agents trained on a rolling window
        lengthHistory = 20;
        
        %DATA GENERATION
        noOverwriteData = false;
        
        %data set sizes
        nTrainingSamples = 20;
        trainingBatch = []; %batches of training data are indexed. This is a list of which batch should go to which agent; useful if you want some agents to have the same training data, otherwise they will be given separate batches.
        
        %informative sensors
        %first row for class 1 and second for class 2
        means = [-3 -2 -1 3 2; 3 2 1 -3 -2];
        deviations = [1 2 1 3 1; 1 2 1 1 3];
        
        cal = ones(1,5);
        bias = zeros(1,5);

        p_c1 = 0.5;

        %non-informative sensors
        noninfMean = 0;
        noninfDev = 5;

        %sensors
        nInfSensors = 5;
        nNoninfSensors = 11;
        
        subs = [];
        
        pMode = 0; %probability of choosing the secondary 'error' sensor mode versus the main 'correct' sensor value mode
        class1ErrorOnly = false;
        
        %probability that a sensor will stop being useful
        pSwitch = 0.015; 
        definiteSwitch = [];
        
        %probability missing labels per agent
        pMissing_a = [];
        
        %default value in case pMissing_a is not large enough for the
        %number of agents
        pMissing = 0; %suggestion: 0.01;
        
        %gaussian-distributed number missing labels
        mean_ml_a = [];
        meanMissingLength = 10; %default
        
        dev_ml_a = [];
        devMissingLength = 2; % default
        
        %probability of corrupt labels for each agent - 
        %corruption of labels is done per-agent, correct label is still
        %recorded
        pCorrupt_a = [];
        pCorrupt = 0; %suggestion: 0.01;
        
        %probability of flipping the labels that an agent sees
        %agent with flipped labels ends up making 100% incorrect
        %predictions
        pFlip_a = [];
        pFlipLabel = 0; % suggestion for moderate flipping: 0.01;
        
        %probability of flipping a sensor - remains informative but agent
        %has to relearn weights for it
        pFlip_s = [];      
        pFlipSensor = 0; %suggestion: 0.01;
        
        %record which labels were corrupted (if any)
        flippedLabels;
        corruptLabels;
        missingLabels;
        
        %probability that a sensor will become "noisy" and start making
        %errors with probability pMode. This is different to "switching"
        %above, where sensors just swap around.
        defNoiseChanges = [];
        switchPerAgent = 0;
        switchDuration = 1;

        %SENSOR ALLOCATION
        reuseSubs = true;
        
        maxSensors = 8;
        minSensors = 3;
    end
    
    methods
        function obj = SynthSettings(topSaveDir)
            if exist('topSaveDir','var')
                obj.topSaveDir = topSaveDir;
            end
        end
        
        function dirName = getDataDir(obj)
            dirName = sprintf('%s/baseData/', obj.topSaveDir);
        end
        
        function nSensors = nSensors(obj)
                nSensors = obj.nInfSensors + obj.nNoninfSensors;
        end        
    end
    
end

