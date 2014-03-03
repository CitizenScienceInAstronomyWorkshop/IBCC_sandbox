function runGZSN(rootDir)

display('starting GZSN classifier combination');

reusePriorRun = false;

%If the paths are added to the default MATLAB path we can comment these out
%and save a bit of time
if nargin<1
    rootDir = '/homes/49/edwin/matlab/combination/gzsnBundle/';
end
addpath([rootDir '../../../matlab']);
dataDirectory = [rootDir 'data'];

%filename for the dataset should be here. Requires the following fields:
%1. classification ID
%2. agent ID
%3. asset ID
%4. PTF type
%5. PTF class
filename = [dataDirectory '/gzsnData.csv'];

%output directory
allFoldsDir = [dataDirectory '/gzsnCombinedScores/'];  
saveAsCsv = true;

%filter out data with not enough responses from frequently-responding agents
minAgentResp = 50; %min no. responses for a frequently-responding agent 
minFreqAgents = 10; %min number of freqently-responding agents per data point
startNegAssets = 1;
maxNoAssets = 50000;

if ~reusePriorRun || ~exist('pointsToDrop','var')
    pointsToDrop = []; %Points we have already processed. The results from these points are in the priors.
else
    pointsToDrop = nonZeroLabels;
    clear idxsToKeep;
end
%Set to false when running live, to true when testing. Save repeating the
%pre-processing step if the raw input data has not changed.
keepDataBetweenRuns = true;

%Set to false to prevent reloading and reprocessing data.
preProcData = true;

%Load the raw dataset - only activated if preProcData is also true.
loadDataFromFile = true;

%Set to false to use existing division of data into n folds; set to true to
%repartition data.
recreateFolds = false;

%Treat the labels as data from a reliable agent.
labelsAsAgent = false;

%include unlabelled/untestable data points
includeUnlabelled = false;

nFolds = 1;

%Settings for graphing the error rates after each fold.
drawGraphs = false;
sortResults = false;

%pick the combination methods to use.
combMethods = {...
    combiners.MeanDecision.shortLabel,...  
...%     combiners.SimpleMajorityVoting.shortLabel, ...
...%     combiners.weighted.WeightedMajority.subShortLabel, ...
...%     combiners.bcc.IbccVb.shortLabel,...   
    combiners.bcc.DynIbccVb.shortLabelInitAll,...
...%     combiners.bcc.IbccVb.shortLabelSeq, ...      
...%     combiners.bcc.Ibcc.shortLabel, ...
...%     combiners.weighted.WeightedSum.shortLabel,...
    };

settings.gz.gz_sn1;

if preProcData && ~exist('snBaseOutputs', 'var')
    if ~exist('snRawData', 'var') 
        snRawData = [];
    end
    [snBaseOutputs, snRawData, labels, typeLabels, typeAssets, assetIds] = ...
        reloadGZSNData(true, true, filename, snRawData);
    
    nTypeKnown = length(typeLabels);
    nAgents = max(snBaseOutputs{1});
    nAssets = max(snBaseOutputs{2});
end

if reusePriorRun && exist('agentRatings','var')
    Alpha = agentRatings{1}{3};
else
    Alpha = [];
end

% Subsample the data to rebalance class proportions. Assuming positive 
% examples are fewer than negative ones, get rid of some negative examples.
if isempty(pointsToDrop)
    if ~exist('idxsToKeep', 'var')
%         [idxsToKeep, subLabels] = dataPrep.balanceClasses(snBaseOutputs, ...
%             labels, minAgentResp, minFreqAgents, includeUnlabelled, true, maxNoAssets);

          %if we don't want any subsampling, use this
          subLabels = labels;
          idxsToKeep = ismember(snBaseOutputs{2}, find(labels~=0));
    end
    testData = {};
    testData{1} = snBaseOutputs{1}(idxsToKeep);
    testData{2} = snBaseOutputs{2}(idxsToKeep);
    testData{3} = snBaseOutputs{3}(idxsToKeep);
else
    idxsToDrop = ismember(snBaseOutputs{2}, pointsToDrop);
        
    testData{1} = snBaseOutputs{1}(~idxsToDrop);
    testData{2} = snBaseOutputs{2}(~idxsToDrop);
    testData{3} = snBaseOutputs{3}(~idxsToDrop);    
    
    assetsToFilter = setDiff(find(labels~=0), pointsToDrop);
    agentsAlreadySeen = unique(snBaseOutputs{1}(idxsToDrop));
    agentValidity = sparse(max(testData{1}), 1);
    agentValidity(~ismember(testData{1}, agentsAlreadySeen)) = 1;
    
    [idxsToKeep, subLabels] = dataPrep.filterData(testData, assetsToFilter, ...
        [], agentValidity, labels, 1, 1, true);
    
    testData{1} = snBaseOutputs{1}(idxsToKeep);
    testData{2} = snBaseOutputs{2}(idxsToKeep);
    testData{3} = snBaseOutputs{3}(idxsToKeep);    
end

nonZeroLabels = find(subLabels~=0);

if (nFolds > 1) && (recreateFolds || ~exist('c_class', 'var'))
    %run with k-fold cross validation
    nLabels = numel(nonZeroLabels);
    c_class = cvpartition(nLabels, 'kfold', nFolds);
end

%If the labels will be supplied as data, create a new agent ID for the
%label provider. Two sets of labels: PTF type and PTF class.
if labelsAsAgent
    ptfTypeAgent = max(snBaseOutputs{1}) + 1;
    ptfClassAgent = ptfTypeAgent + 1;
end

%record the results of iterative methods at each iteration
progressMonitor = scoring.CombinerProgressMonitor([], length(combMethods));
progressMonitor.method = 'absolute error';

totalTime = zeros(1, length(combMethods));    

resultsAllFolds = []; %results of combination for all n folds
labelsAllFolds = []; %concatenate labels for all n folds
resultTable = []; %results of combination with asset IDs in first column

agentRatings = cell(nFolds, 1);

for i=1:nFolds

    display(['Fold: ' num2str(i)]);
    
    foldData = cell(length(testData), 1);

    if nFolds > 1
        testIdxs = nonZeroLabels(c_class.test(i));
        trainingIdxs = nonZeroLabels(c_class.training(i));      
        typeTrainIdxs = ~ismember(typeAssets, testIdxs);
        typeTrainAssets = typeAssets(typeTrainIdxs);
    else
        uniqueAssets = unique(testData{2});
        testAssets = zeros(size(subLabels,1), size(subLabels,2))-1;
        testAssets(uniqueAssets) = subLabels(uniqueAssets);
        testIdxs = nonZeroLabels;

        display(['number of test data points: ' num2str(length(testIdxs))]);        
       
        trainingIdxs = nonZeroLabels;
        typeTrainIdxs = ~ismember(typeAssets, testIdxs);
        typeTrainAssets = typeAssets(typeTrainIdxs);
    end

    %data 1 is agents
    if labelsAsAgent   
        foldData{1} = [testData{1}; ... %data from human classifiers
                ptfTypeAgent*ones(length(typeTrainIdxs), 1); ... %type data from ptf (assets with a ptf class)
            ];%ptfClassAgent*ones(sum(c_class.training(i)), 1)]; %class data from ptf

        foldData{2} = [testData{2}; ...%data point no.
                typeTrainAssets];%classAssets(c_class.training(i))];

        foldData{3} = [testData{3}; ...%scores given by agents
                typeLabels(typeTrainIdxs)];%classLabels(c_class.training(i))];
        
        trainingLabels = zeros(size(subLabels,1),size(subLabels,2));
    else
        foldData = testData;
        trainingLabels = subLabels;
        %if we are using all indices as test indices, also use them all as
        %training indices.
        if nFolds>1   
            trainingLabels(testIdxs) = 0;
        end
    end

    testLabels = subLabels(testIdxs);
    
    display(['length of fold data=' num2str(length(foldData{1}))]);
    runner = GzRunner(expSettings, foldData, trainingLabels, subLabels, Alpha);
    
    if nFolds > 1 %if nFolds==1 we are not splitting into folds
        progressMonitor.updateLabels(subLabels, testIdxs, length(trainingIdxs));
        runner.progressMonitor = progressMonitor;
    end

    [results, agentRatings{i}] = runner.runCombineAll(false, 'dontTouch', drawGraphs, combMethods, sortResults); 

    totalTime = totalTime + runner.runTime;

    testResults = results(:, testIdxs);        
    
    resultsAllFolds = cat(2, resultsAllFolds, testResults);
    labelsAllFolds = cat(2, labelsAllFolds, testLabels');

    if ~exist(sprintf('%sall%dFolds', allFoldsDir, nFolds), 'dir');
        mkdir(sprintf('%sall%dFolds', allFoldsDir, nFolds));
    end

    %save results after each fold
    if saveAsCsv
        %map test Idxs back to assetIds
        testAssetIds = assetIds(testIdxs);
        resultTable = cat(1, resultTable, [testAssetIds testResults']);
        dlmwrite(sprintf('%s%s_minPosR%d_minFreqA%d.csv', allFoldsDir, ...
            datestr(now, 'yy_mm_dd__HH_MM_SS'), minAgentResp, minFreqAgents), ...
            resultTable);
    else
        save(sprintf(...
            '%sall%dFolds/combin_minPosR%d_minFreqA%d.mat', ...
            allFoldsDir, nFolds, minAgentResp, minFreqAgents), 'resultsAllFolds');    
        save(sprintf(...
            '%sall%dFolds/labels_minPosR%d_minFreqA%d.mat', ...
            allFoldsDir, nFolds, minAgentResp, minFreqAgents), 'labelsAllFolds');
        save(sprintf(...
            '%sall%dFolds/cvClass_minPosR%d_minFreqA%d.mat', ...
            allFoldsDir, nFolds, minAgentResp, minFreqAgents), 'c_class');      
        save(sprintf(...
            '%sall%dFolds/progMon_minPosR%d_minFreqA%d.mat', ...
            allFoldsDir, nFolds, minAgentResp, minFreqAgents), 'progressMonitor');         
    end
end

avgTime = totalTime ./ nFolds;
display(['Average time per fold for each combination method: ' num2str(avgTime)]);

if ~keepDataBetweenRuns
    clear snRawData
    clear snBaseOutputs
    clear idxsToKeep
end
display('completed classifier combination.');
end
