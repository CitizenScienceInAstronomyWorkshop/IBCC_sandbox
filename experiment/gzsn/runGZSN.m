%move this line out to make this not specific to GZSN! Put at the top of
%the experiment-specific run scripts. This way we avoid overwriting any 
%experiment-specific settings using this script.
% settings.gz.gz_sn1;

%Runs tests on a dataset using n-fold cross validation or running over the
%whole dataset using all labels provided.

%SETTINGS -----------------------------------------------------------------

%Default settings objects are created unless these objects already exist.
%Provide experiment-specific settings by setting properties of the 
%following objects before calling this script
% ExpSettings: expSet
% BccSettings: bccSet
% cell array of strings: combMethods

if ~exist('expSet', 'var')
    expSet = settings.ExpSettings();
end

display(['root dir for this application: ' expSet.rootDir]);
display(['data dir: ' expSet.dataDir]);
display(['output dir: ' expSet.outputDir]);

if ~exist('bccSet', 'var')
    bccSet = settings.BccSettings();
end

%pick the combination methods to use.
if ~exist('combMethods','var')
    combMethods = {...
        combiners.MeanDecision.shortLabel,...  
        combiners.bcc.IbccVb.shortLabel,...       
        };
end

%--------------------------------------------------------------------------

display('starting GZSN classifier combination');

if ~expSet.reusePriorRun || ~exist('pointsToDrop','var')
    pointsToDrop = []; %Points we have already processed. The results from these points are in the priors.
else
    pointsToDrop = nonZeroLabels;
    clear idxsToKeep;
end

if expSet.preProcData && ~exist('snBaseOutputs', 'var')
    if ~exist('snRawData', 'var') 
        snRawData = [];
    end
    [snBaseOutputs, snRawData, labels, typeLabels, typeAssets, assetIds] = ...
        reloadGZSNData(true, true, expSet.inputFile, snRawData);
    
    nTypeKnown = length(typeLabels);
    nAgents = max(snBaseOutputs{1});
    nAssets = max(snBaseOutputs{2});
end

if expSet.reusePriorRun && exist('agentRatings','var')
    Alpha = agentRatings{1}{3};
% elseif ~exist('Alpha', 'var') %keep alpha otherwise
else
    Alpha = bccSet.Alpha;
% else
    display('using pre-set Alpha hyperparameters');
end

% Subsample the data to rebalance class proportions. Assuming positive 
% examples are fewer than negative ones, get rid of some negative examples.
if isempty(pointsToDrop)
    if ~exist('idxsToKeep', 'var')
        [idxsToKeep, subLabels] = dataPrep.balanceClasses(snBaseOutputs, ...
            labels, expSet.minAgentResp, expSet.minFreqAgents, ...
            expSet.includeUnlabelled, true, expSet.maxNoAssets);
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

if (expSet.recreateFolds && expSet.nFolds > 1) || ~exist('c_class', 'var')
    %run with k-fold cross validation
    nLabels = numel(nonZeroLabels);
    c_class = cvpartition(nLabels, 'kfold', expSet.nFolds);
end

%If the labels will be supplied as data, create a new agent ID for the
%label provider. Two sets of labels: PTF type and PTF class.
if bccSet.typeLabelsAsAgent
    ptfTypeAgent = max(snBaseOutputs{1}) + 1;
%     ptfClassAgent = ptfTypeAgent + 1;
end

%record the results of iterative methods at each iteration
progressMonitor = scoring.CombinerProgressMonitor([], length(combMethods));
progressMonitor.method = 'absolute error';

totalTime = zeros(1, length(combMethods));    

resultsAllFolds = []; %results of combination for all n folds
labelsAllFolds = []; %concatenate labels for all n folds
resultTable = []; %results of combination with asset IDs in first column

agentRatings = cell(expSet.nFolds, 1);

for i=1:expSet.nFolds

    display(['Fold: ' num2str(i)]);
    
    foldData = cell(length(testData), 1);

    if expSet.nFolds > 1
        testIdxs = nonZeroLabels(c_class.test(i));
        trainingIdxs = nonZeroLabels(c_class.training(i));     
        
        typeTrainIdxs = ~ismember(typeAssets, testIdxs);
        typeTrainAssets = typeAssets(typeTrainIdxs);
    else
        uniqueAssets = unique(testData{2});
        testAssets = zeros(size(subLabels,1), size(subLabels,2))-1;
        testAssets(uniqueAssets) = subLabels(uniqueAssets);
        testIdxs = find(testAssets==0);
       
        trainingIdxs = nonZeroLabels;
        typeTrainIdxs = ~ismember(typeAssets, testIdxs);
        typeTrainAssets = typeAssets(typeTrainIdxs);
    end

    %data 1 is agents
    if bccSet.typeLabelsAsAgent && ~isempty(typeTrainIdxs)
        foldData{1} = [testData{1}; ... %data from human classifiers
                ptfTypeAgent*ones(sum(typeTrainIdxs), 1); ... %type data from ptf (assets with a ptf class)
            ];%ptfClassAgent*ones(sum(c_class.training(i)), 1)]; %class data from ptf

        foldData{2} = [testData{2}; ...%data point no.
                typeTrainAssets];%classAssets(c_class.training(i))];

        foldData{3} = [testData{3}; ...%scores given by agents
                typeLabels(typeTrainIdxs)];%classLabels(c_class.training(i))];
        
        trainingLabels = zeros(size(subLabels,1),size(subLabels,2));
        bccSet.targetsAsAgents = 1;
    else
        foldData = testData;
        trainingLabels = subLabels;
        trainingLabels(testIdxs) = 0;
        bccSet.targetsAsAgents = [];
        bccSet.trustedAlpha = [];       
    end

    testLabels = subLabels(testIdxs);
    
    display(['length of fold data=' num2str(length(foldData{1}))]);
    runner = GzRunner(expSet, foldData, trainingLabels, subLabels, ...
        Alpha, bccSet.targetsAsAgents, bccSet.trustedAlpha);
    
    if nFolds > 1 %if nFolds==1 we are not splitting into folds
        progressMonitor.updateLabels(subLabels, testIdxs, length(trainingIdxs));
        runner.progressMonitor = progressMonitor;
    end

    [results, agentRatings{i}] = runner.runCombineAll(false, 'dontTouch', drawGraphs, combMethods, sortResults); 

    graphCombinationInProgress( testIdxs, results(:,testIdxs), testLabels, foldData);
    
    totalTime = totalTime + runner.runTime;

    testResults = results(:, testIdxs);        
    
    resultsAllFolds = cat(2, resultsAllFolds, testResults);
    labelsAllFolds = cat(2, labelsAllFolds, testLabels');

    if ~exist(sprintf('%s/all%dFolds', expSet.outputDir, expSet.nFolds), 'dir');
        mkdir(sprintf('%s/all%dFolds', expSet.outputDir, expSet.nFolds));
    end

    %save results after each fold
    if expSet.saveAsCsv
        
        %map test Idxs back to assetIds
        testAssetIds = assetIds(testIdxs);
        resultTable = cat(1, resultTable, [testAssetIds testResults']);
        dlmwrite(sprintf('%s/%s_minPosR%d_minFreqA%d.csv', expSet.outputDir, ...
            datestr(now, 'yy_mm_dd__HH_MM_SS'), expSet.minAgentResp, expSet.minFreqAgents), ...
            resultTable);
    else
        save(sprintf(...
            '%sall%dFolds/combin_minPosR%d_minFreqA%d.mat', ...
            expSet.outputDir, expSet.nFolds, expSet.minAgentResp, expSet.minFreqAgents), 'resultsAllFolds');    
        save(sprintf(...
            '%sall%dFolds/labels_minPosR%d_minFreqA%d.mat', ...
            expSet.outputDir, expSet.nFolds, expSet.minAgentResp, expSet.minFreqAgents), 'labelsAllFolds');
        save(sprintf(...
            '%sall%dFolds/cvClass_minPosR%d_minFreqA%d.mat', ...
            expSet.outputDir, expSet.nFolds, expSet.minAgentResp, expSet.minFreqAgents), 'c_class');      
        save(sprintf(...
            '%sall%dFolds/progMon_minPosR%d_minFreqA%d.mat', ...
            expSet.outputDir, expSet.nFolds, expSet.minAgentResp, expSet.minFreqAgents), 'progressMonitor');         
    end
end

avgTime = totalTime ./ expSet.nFolds;
display(['Average time per fold for each combination method: ' num2str(avgTime)]);

if ~expSet.keepDataBetweenRuns
    clear snRawData
    clear snBaseOutputs
    clear idxsToKeep
end
display('completed classifier combination.');