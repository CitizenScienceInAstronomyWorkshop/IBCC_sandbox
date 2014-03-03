classdef ExpSettings < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        %CHARACTERISTICS OF APPLICATION------------------------------------
        
        nScores = 2; %possible base classifier outputs
        nClasses = 2; %possible classes
        
        %value in results when a base classifier hasn't scored a sample
%         noScore = -1;        
        noScore = 0;        
        
        nSamples = 0; %number of data points to classify
        
        spreadTargets = true; %when choosing known targets randomly, 
        %spread them out or take first few
                
        %GRAPH OUTPUT------------------------------------------------------
        
        %Settings for graphing the error rates after each fold.
        drawGraphs = false;
        sortResults = false;
        saveImgs = true;
        graphVisible = 'On';
        
        %DATA DIR/FILENAMES------------------------------------------------
        topSaveDir = 'unknownSaveDirectory';
        expLabel = 'unknownExperiment';
        multCombTestFile = 'multCombinerTests.mat';
                        
        nRepeats = 5;
        
        %DATA INPUT AND OUTPUT---------------------------------------------
        rootDir = '/homes/49/edwin/matlab/distrib/gzsnBundle/';    
        dataDir;
        
        %filename for the dataset should be here. Requires the following fields:
        %1. classification ID
        %2. agent ID
        %3. asset ID
        %4. PTF type
        %5. PTF class
        %6. score
        inputFile;
        
        outputDir;
        saveAs = 'csv';
        
        %write over data from previous runs of the experiment
        dontOverwrite = false;
        
        %include the training data points in the results. Useful if writing 
        %final output e.g. galaxy zoo needs an entry for all assets.        
        writeTrainingResults = false;
        
        %SUBSAMPLING...----------------------------------------------------
        %data points by specifying minimal requirements for base
        %classifiers 
        
        %filter out data with not enough responses from frequently-responding agents
        minAgentResp = 1; %min no. responses for a frequently-responding agent 
        minFreqAgents = 1; %min number of freqently-responding agents per data point
        startNegAssets = 1;
        maxNoAssets = 1000000000;       
        
        %include unlabelled/untestable data points
        includeUnlabelled = false;
        
        %VOTING METHODS ---------------------------------------------------
        
        voteThreshold = 0.5;
        postScores = true;        
        
        %N-FOLD CROSS VALIDATION-------------------------------------------
        
        %Set to false to use existing division of data into n folds; 
        %set to true to repartition data.
        recreateFolds = false;
        
        nFolds = 5; 
        
        %instead of folds, split the unlabelled data randomly into batches 
        %and use all training data available - useful if algorithm 
        %performance drops with many unlabelled data points
        batches = false; 
        
        %ALTERNATIVE TO CROSS VALIDATION ----------------------------------
        %For when we have very few labels, we might not want to run every
        %fold. Also, we want to run same number of times for each
        %proportion of known labels
        propKnown = [0.1 0.2 0.3];
        iPropKnown = 1;
        
        %USE DATA FROM PREVIOUS RUNS?--------------------------------------
        
        %resuse previous results to set priors for this run
        reusePriorRun = false;
        
        %Delete 
        %Set to false when running live, to true when testing. Save 
        %repeating the pre-processing step if the raw input data has not 
        %changed.
        keepDataBetweenRuns = true;

        %Set to false to prevent reloading and reprocessing data.
        preProcData = true;

        %Load the raw dataset - only activated if preProcData is also true.
        loadDataFromFile = true;
    end
    
    properties (Constant)
        dlrPath = '~/svn/edwin/matlab/third_party/dlr';
    end
    
    methods
        function obj = ExpSettings(rootDir, inputFilename, outputSubdir)
            %each cell contains posteriors: each set of results is a row,
            %each row is a repeat with same settings.
            
            if ~exist('inputFilename','var')
                inputFilename = '/baseData.csv';
            end
            if ~exist('outputSubdir','var')
                outputSubdir = '/output';
            end
            
            obj.setDirNames(rootDir, inputFilename, outputSubdir);
        end
        
        
        function setDirNames(obj, rootDir, inputFilename, outputSubdir)

            if nargin >= 1 && ~isempty(rootDir)
                obj.rootDir = rootDir;
                obj.topSaveDir = rootDir;
            end
            
            if ~exist('inputFilename','var')
                inputFilename = '/baseData.csv';
            end
            if ~exist('outputSubdir','var')
                outputSubdir = '/output';
            end
            
            addpath([obj.rootDir '../../']);
            obj.dataDir = [obj.rootDir '/data/'];
            obj.inputFile = [obj.dataDir inputFilename];
            obj.outputDir =  [obj.rootDir '/' outputSubdir];
            
            if ~exist(obj.outputDir, 'dir')
                mkdir(obj.outputDir);
            end
        end
        
        function dirName = getDataDir(obj)
            dirName = obj.dataDir;
            if ~exist(dirName,'dir')
                mkdir(dirName);
            end
        end
        
        function dirName = getCombinerDir(obj)
            dirName = sprintf('%s/%s/', obj.topSaveDir, obj.expLabel);
            if ~exist('dirName','dir')
               mkdir(dirName);
            end
        end      
    end
    
end

