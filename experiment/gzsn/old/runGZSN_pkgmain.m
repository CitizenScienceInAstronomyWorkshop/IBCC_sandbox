function [results labels agentRatings ptfTypeAgent] = runGZSN_pkgmain(rootDir)

display('starting GZSN classifier combination');

settings.gz_sn1;

%If the paths are added to the default MATLAB path we can comment these out
%and save a bit of time
if nargin<1
    expSet.rootDir = '/homes/49/edwin/matlab/combination/gzsnBundle/';
    addpath([rootDir '../../../matlab']);
end
dataDirectory = [rootDir 'data'];

%filename for the dataset should be here. Requires the following fields:
%1. classification ID
%2. agent ID
%3. asset ID
%4. PTF type
%5. PTF class
filename = [dataDirectory '/input.csv'];

%output directory
allFoldsDir = [dataDirectory '/gzsnCombinedScores/'];  

%Treat the *type* labels as data from a reliable agent.
labelsAsAgent = true;

%Settings for graphing the error rates after each fold.
drawGraphs = false;
sortResults = false;

%pick the combination methods to use.
combMethods = {...
...%    combiners.MeanDecision.shortLabel,...  
...%     combiners.SimpleMajorityVoting.shortLabel, ...
...%     combiners.weighted.WeightedMajority.subShortLabel, ...
     combiners.bcc.IbccVb.shortLabel,...   
    };


%turns off the assumption that no PTF type --> not supernova
[testData, snRawData, labels, typeLabels, typeAssets, assetIds] = ...
    reloadGZSNData(true, true, filename, [], false);

nTypeKnown = length(typeLabels)
nAgents = max(testData{1})
nAssets = max(testData{2})

nonZeroLabels = find(labels~=0);

%If the labels will be supplied as data, create a new agent ID for the
%label provider. Two sets of labels: PTF type and PTF class.
if bccSet.screenLabelsAsAgent
    ptfTypeAgent = max(testData{1}) + 1;
end

%record the results of iterative methods at each iteration
progressMonitor = scoring.CombinerProgressMonitor([], length(combMethods));
progressMonitor.method = 'absolute error';

totalTime = zeros(1, length(combMethods));    

resultTable = []; %results of combination with asset IDs in first column

foldData = cell(length(testData), 1);

display(['number of unknown data points: ' num2str(nAssets-length(nonZeroLabels))]);        

%data 1 is agents
if bccSet.screenLabelsAsAgent && ~isempty(typeAssets)
    foldData{1} = [testData{1}; ... %data from human classifiers
            ptfTypeAgent*ones(length(typeAssets), 1); ... %type data from ptf (assets with a ptf class)
        ];

    foldData{2} = [testData{2}; ...%data point no.
            typeAssets];

    foldData{3} = [testData{3}; ...%scores given by agents
            typeLabels];

    trainingLabels = labels;
    targetsAsAgents = 1;
    trustedAlpha = [50 10 1; 1 20 40];
else
    foldData = testData;
    trainingLabels = labels;
    targetsAsAgents = [];
    trustedAlpha = [];
end

runner = ExpRunner(expSet, bccSet, [], [], labels, foldData, trainingLabels);
[results, agentRatings] = runner.runCombineAll(false, 'dontTouch', drawGraphs, combMethods, sortResults); 

totalTime = totalTime + runner.runTime;     

if ~exist(sprintf('%s', allFoldsDir), 'dir');
    mkdir(sprintf('%s', allFoldsDir));
end

%map test Idxs back to assetIds
resultTable = cat(1, resultTable, [assetIds results']);
outFile = sprintf('%sibccOutputs.csv', allFoldsDir);
display(['writing to ' outFile]);
dlmwrite(outFile, resultTable);

avgTime = totalTime;
display(['time for each combination method: ' num2str(avgTime)]);

display('completed classifier combination.');
end
