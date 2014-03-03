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

if expSet.reusePriorRun && exist('agentRatings','var')
    Alpha = agentRatings{1}{3};
else
    Alpha = bccSet.Alpha;
    display('using pre-set Alpha hyperparameters');
end

%LOADING DATA -------------------------------------------------------------
%...if not already done - avoid repeating preprocessing large files
if expSet.preProcData && ~exist('baseOutputs', 'var')
    if ~exist('rawData', 'var') 
        rawData = [];
    end
    
    if ~exist('sepPosLabFile','var')
        sepPosLabFile = [];
    end
    
    if ~exist('sepNegLabFile','var')
        sepNegLabFile = [];
    end  
    
    if ~exist('rawColumnMap', 'var')
        rawColumnMap = [];
    end
    
    if ~exist('columnString', 'var')
        columnString = [];
    end
    
    if ~exist('legalScores', 'var')
        legalScores = [1 2 3];
    end
    
    if ~exist('interpretEmptyScreenLabels', 'var')
        interpretEmptyScreenLabels = true;
    end
    
    [baseOutputs, rawData, labels, screenLabels, screenObjs, objIds, agentIds] = ...
        loadBaseData(expSet.inputFile, rawData, sepPosLabFile, ...
            sepNegLabFile, rawColumnMap, columnString, legalScores, interpretEmptyScreenLabels);
    
    nTypeKnown = length(screenLabels);
    nAgents = max(baseOutputs{1});
    nAssets = max(baseOutputs{2});
end

%If the labels will be supplied as data, create a new agent ID for the
%label provider. Two sets of labels: PTF type and PTF class.
if bccSet.screenLabelsAsAgent
    screenAgent = max(baseOutputs{1}) + 1;
end

%SUBSAMPLING -------------------------------------------------------------- 

if ~expSet.reusePriorRun || ~exist('pointsToDrop','var')
    %Points we have already processed, in this case none - process all
    pointsToDrop = []; 
else
    %The results from these points are in the priors - don't process again       
    pointsToDrop = nonZeroLabels;
    %Have to clear this as we are working with a new subsample. Normally, 
    %retain idxsToKeep between runs to avoid repeating subsampling step.
    clear idxsToKeep;
end

% Subsample the data to rebalance class proportions. Assuming positive 
% examples are fewer than negative ones, get rid of some negative examples.
if isempty(pointsToDrop)
    
    if ~exist('keepInfreqAgents','var')
        keepInfreqAgents = true;
    end
    
    if ~exist('filterPos','var')
        filterPos = false;
    end
    
    if ~exist('idxsToKeep', 'var')
        [idxsToKeep, subLabels] = dataPrep.balanceClasses(baseOutputs, ...
            labels, expSet.minAgentResp, expSet.minFreqAgents, ...
            expSet.includeUnlabelled, keepInfreqAgents, expSet.maxNoAssets, ...
            filterPos);%, expSet.startNegAssets);
        
%         nSamples = floor(length(baseOutputs{1})/4);

%         idxsToKeep = (expSet.startNegAssets*nSamples - nSamples+1):expSet.startNegAssets*nSamples;
%         subLabels = sparse(length(labels), 1);
%         assetsToKeep = baseOutputs{2}(idxsToKeep);
%         subLabels(assetsToKeep) = labels(assetsToKeep);
    end
    testData = {};
    testData{1} = baseOutputs{1}(idxsToKeep);
    testData{2} = baseOutputs{2}(idxsToKeep);
    testData{3} = baseOutputs{3}(idxsToKeep);
else
    idxsToDrop = ismember(baseOutputs{2}, pointsToDrop);
        
    testData{1} = baseOutputs{1}(~idxsToDrop);
    testData{2} = baseOutputs{2}(~idxsToDrop);
    testData{3} = baseOutputs{3}(~idxsToDrop);    
    
    assetsToFilter = setDiff(find(labels~=0), pointsToDrop);
    agentsAlreadySeen = unique(baseOutputs{1}(idxsToDrop));
    agentValidity = sparse(max(testData{1}), 1);
    agentValidity(~ismember(testData{1}, agentsAlreadySeen)) = 1;
    
    [idxsToKeep, subLabels] = dataPrep.filterData(testData, assetsToFilter, ...
        [], agentValidity, labels, 1, 1, true);
    
    testData{1} = baseOutputs{1}(idxsToKeep);
    testData{2} = baseOutputs{2}(idxsToKeep);
    testData{3} = baseOutputs{3}(idxsToKeep);    
end

nonZeroLabels = find(subLabels~=0);
uniqueAssets = unique(testData{2});
zeroLabels = find(subLabels==0);
zeroLabels = zeroLabels(ismember(zeroLabels, uniqueAssets));

display(['Complete dataset contains ' num2str(length(nonZeroLabels)) ' known labels and ' num2str(length(zeroLabels)) ' unknown labels. Total = ' num2str(length(subLabels))]);
display(['There are ' num2str(length(unique(testData{1}))) ' agents']);
resultsAllFolds = []; %results of combination for all n folds
labelsAllFolds = []; %concatenate labels for all n folds
testIdxsAllFolds = []; %indexes into the results for non-training samples (only makes sense when we include training samples in the results)
objIdsAllFolds = [];
resultTable = []; %results of combination with asset IDs in first column

%PROGRESS MONITOR ---------------------------------------------------------

%record the results of iterative methods at each iteration
if ~exist('progressMonitor','var')
    progressMonitor = scoring.CombinerProgressMonitor([], length(combMethods));
    progressMonitor.method = 'AUC';
end

totalTime = zeros(1, length(combMethods));    

%PARTITIONING FOR CROSS VALIDATION ----------------------------------------

if (expSet.recreateFolds || ~exist('c_class', 'var') || c_class.NumTestSets~=expSet.nFolds) && expSet.nFolds > 1
    %run with k-fold cross validation
    nLabels = numel(nonZeroLabels);
    c_class = cvpartition(nLabels, 'kfold', expSet.nFolds);
    
    nRuns = expSet.nFolds;
elseif expSet.nFolds <=1
    nRuns = 1;
end
agentRatings = cell(expSet.nFolds, 1);

%PARTITIONING THE UNLABELLED DATA -----------------------------------------
if expSet.batches 
    if expSet.nFolds>1
        nBatches = expSet.nFolds;
        batchSize = ceil(length(zeroLabels)./expSet.nFolds);
    else
        batchSize = round(length(nonZeroLabels ./ 2));
        nBatches = ceil(length(zeroLabels)./batchSize);
        nRuns = nBatches;
    end
else
    nRuns = expSet.nFolds;
    if nRuns<1
        nRuns = 1;
    end
    nBatches = 1;
    batchSize = length(zeroLabels);
end

batches = cell(1,nBatches);
prevEnd = 0;

for i=1:nBatches
    currentEnd = prevEnd+batchSize;
    if currentEnd > length(zeroLabels)
        currentEnd = length(zeroLabels);
    end
    batches{i} = zeroLabels(prevEnd+1:currentEnd);
    prevEnd = currentEnd;
end

%RUN N-FOLD CROSS VALIDATION ----------------------------------------------

display('starting classifier combination');

aucInt = cell(1, nRuns);
aucIntAlt = cell(1, nRuns);
for i=1:nRuns

    display(['Fold: ' num2str(i)]);
    
    foldData = cell(length(testData), 1);
    %TRAINING/TEST INDEXES FOR THIS PARTITION -----------------------------
    if expSet.nFolds > 1
        if nBatches > 1
            testIdxs = [nonZeroLabels(c_class.test(i)); batches{i}];
        else
            testIdxs = [nonZeroLabels(c_class.test(i)); batches{1}]; 
        end
        trainingIdxs = nonZeroLabels(c_class.training(i));     
    elseif expSet.nFolds==0
        testIdxs = [nonZeroLabels; batches{i}];
        trainingIdxs = [];        
    else %real run in a live application - just want the results, no evaluation
%         uniqueAssets = unique(testData{2});
%         testAssets = zeros(size(subLabels,1), size(subLabels,2))-1;
%         testAssets(uniqueAssets) = subLabels(uniqueAssets);
        testIdxs = batches{i};
        trainingIdxs = nonZeroLabels;
    end
    
    %SCREEN LABELS---------------------------------------------------------
    %E.g. the PTF type in GZSN. Assume these labels are another base classifier. 
    % Commented out code is the same idea for true label class; however, we
    %usually believe class labels are equivalent to the ground truth, so 
    %don't treat them as unreliable agents.
    
    runDataIdxs = ismember(testData{2}, [trainingIdxs; testIdxs]);
    typeTrainIdxs = ~ismember(screenObjs, testIdxs);            
    if bccSet.screenLabelsAsAgent && ~isempty(typeTrainIdxs)
        typeTrainIdxs = typeTrainIdxs(1);
        typeTrainAssets = screenObjs(typeTrainIdxs);            
        
        foldData{1} = [testData{1}(runDataIdxs); ... %data from base classifiers
            screenAgent*ones(sum(typeTrainIdxs), 1);...%labels from screening process
        ];%ptfClassAgent*ones(sum(c_class.training(i)), 1)]; %class labels

        foldData{2} = [testData{2}(runDataIdxs); ...%data point no.
            typeTrainAssets];
            %classAssets(c_class.training(i))];

        foldData{3} = [testData{3}(runDataIdxs); ...%scores given by agents
            screenLabels(typeTrainIdxs)];
            %classLabels(c_class.training(i))];
        
        bccSet.targetsAsAgents = 1;
    else
        foldData = testData;
        bccSet.targetsAsAgents = [];
        bccSet.trustedAlpha = [];       
    end
    
    trainingLabels = subLabels;
    trainingLabels(testIdxs) = 0;    

    %EXPERIMENT RUNNER ----------------------------------------------------
    display(['Quantity of data to process = ' num2str(length(foldData{1}))]);
    
    runner = ExpRunner(expSet, [], bccSet, [],[], subLabels, foldData, [], trainingLabels);
    
    if expSet.nFolds > 1
        progressMonitor.updateLabels(subLabels, testIdxs, length(trainingIdxs));
        runner.progressMonitor = progressMonitor;
    end

    [results, agentRatings{i}] = runner.runCombineAll(false, 'dontTouch', ...
        expSet.drawGraphs, combMethods, expSet.sortResults); 

    totalTime = totalTime + runner.runTime;

    %RESULT OUPUT ---------------------------------------------------------
    
    if expSet.writeTrainingResults
        testResults = results;
        testLabels = subLabels;
        testResultsIdxs = length(labelsAllFolds)+testIdxs';
        testIdxsAllFolds = [testIdxsAllFolds testResultsIdxs];
        objIdsAllFolds = [objIdsAllFolds 1:size(results,2)];
    else
        testResults = results(:, testIdxs);   
        testLabels = subLabels(testIdxs);   
        testResultsIdxs = length(labelsAllFolds)+1:length(labelsAllFolds)+length(testIdxs);
        testIdxsAllFolds = [testIdxsAllFolds testResultsIdxs];
        objIdsAllFolds = [objIdsAllFolds testIdxs'];
    end
    
    if expSet.drawGraphs
        graphCombinationInProgress(testIdxs, results(:,testIdxs), testLabels, foldData);
    end
        
    %show intermediate results
    if expSet.nFolds>1 || expSet.nFolds==0
        aucOnly = true;
        aucInt{i} = graphs.ClassifierPerformanceGraph.drawRoc(results(:, testIdxs), ...
            subLabels(testIdxs)'-1, combMethods, false, aucOnly, true);
        aucInt{i}

        if exist('labelsBoth', 'var')
            aucIntAlt{i} = graphs.ClassifierPerformanceGraph.drawRoc(results(:, testIdxs), ...
            labelsBoth(testIdxs)'-1, combMethods, false, aucOnly, true);
            aucIntAlt{i}
        end;
    end
    
    resultsAllFolds = cat(2, resultsAllFolds, testResults);
    labelsAllFolds = cat(2, labelsAllFolds, testLabels');

    %save results after each fold
    if strcmp(expSet.saveAs, 'csv')
        
        
        %map test Idxs back to objIds
        if expSet.writeTrainingResults
            testObjIds = objIds;
        else
            testObjIds = objIds(testIdxs);
        end
        
%         resultTable = cat(1, resultTable, [testObjIds testResults']);
        for c=1:length(combMethods)
            if expSet.dontOverwrite
                outFile = sprintf('%s/%s_%s.csv', expSet.outputDir, combMethods{c},...
                    datestr(now, 'yy_mm_dd__HH_MM_SS'));    
            else
                outFile = sprintf('%s/%s_outputs.csv', expSet.outputDir, combMethods{c});    
            end
        
            fid = fopen(outFile, 'w');
            linesPrinted = 0;
            for n=1:size(testResults,2)
                if ~ismember(n,baseOutputs{2})
                    continue;
                end
                string = sprintf('%d, %4f\n', testObjIds(n), testResults(c,n));
                fprintf(fid,string);
                linesPrinted = linesPrinted + 1;
            end
            fclose(fid); 
            display(['Lines saved to file: ' num2str(linesPrinted)]);
        end
                
%         dlmwrite(outFile, resultTable);
    elseif strcmp(expSet.saveAs, 'mat')
        if ~exist(sprintf('%s/all%dFolds', expSet.outputDir, expSet.nFolds), 'dir');
            mkdir(sprintf('%s/all%dFolds', expSet.outputDir, expSet.nFolds));
        end        
        
        save(sprintf(...
            '%sall%dFolds/combin_minPosR%d_minFreqA%d.mat', ...
            expSet.outputDir, expSet.nFolds, expSet.minAgentResp, expSet.minFreqAgents), 'resultsAllFolds');    
        save(sprintf(...
            '%sall%dFolds/labels_minPosR%d_minFreqA%d.mat', ...
            expSet.outputDir, expSet.nFolds, expSet.minAgentResp, expSet.minFreqAgents), 'labelsAllFolds');
        if exist('c_class','var')
            save(sprintf(...
                '%sall%dFolds/cvClass_minPosR%d_minFreqA%d.mat', ...
                expSet.outputDir, expSet.nFolds, expSet.minAgentResp, expSet.minFreqAgents), 'c_class');      
        end
        save(sprintf(...
            '%sall%dFolds/progMon_minPosR%d_minFreqA%d.mat', ...
            expSet.outputDir, expSet.nFolds, expSet.minAgentResp, expSet.minFreqAgents), 'progressMonitor');         
    end
end

avgTime = totalTime ./ nRuns;
display(['Average time per fold for each combination method: ' num2str(avgTime)]);

%TIDY UP AFTER EXPERIMENT IS COMPLETE -------------------------------------

if ~expSet.keepDataBetweenRuns
    clear rawData
    clear baseOutputs
    clear idxsToKeep
end
display('completed classifier combination.');

if expSet.nFolds==1
    return
end

%DRAW A LOVELY ROC CURVE --------------------------------------------------
auc = graphs.ClassifierPerformanceGraph.drawRoc(resultsAllFolds(:,testIdxsAllFolds), ...
    labelsAllFolds(testIdxsAllFolds)-1, combMethods, false, false, true)
lam = cell(1,length(combMethods));

Cmat = sparse(baseOutputs{1}, baseOutputs{2}, baseOutputs{3});

avgAuc = zeros(1, length(combMethods));
brier = zeros(1, length(combMethods));

for c=1:length(combMethods)
    %assume class 2 is positive and class 1 is negative
    fp = resultsAllFolds(c,testIdxsAllFolds)>=0.5 & labelsAllFolds(testIdxsAllFolds)==1;
%     resultsAllFolds(c,testIdxsAllFolds(fp))
%     Cmat(:, objIdsAllFolds(testIdxsAllFolds(fp)))
    fn = resultsAllFolds(c,testIdxsAllFolds)<0.5 & labelsAllFolds(testIdxsAllFolds)==2;
%     resultsAllFolds(c,testIdxsAllFolds(fn))
%     Cmat(:, objIdsAllFolds(testIdxsAllFolds(fn)))

    nFp = sum(fp);
    nFn = sum(fn);

    nNeg = sum(labelsAllFolds(testIdxsAllFolds)==1);
    nPos = sum(labelsAllFolds(testIdxsAllFolds)==2);

    if nNeg==0
        fpr = 0.5;
    else   
        fpr = nFp ./ nNeg;
    end

    if nPos==0
        fnr = 0.5;
    else
        fnr = nFn ./ nPos;
    end

    lamc = (log(fpr./(1-fpr)) + log(fnr./(1-fnr)))./2;
    lam{c} = exp(lamc) ./ (1+exp(lamc));

    display([combMethods{c} ': LAM = ' num2str(lam{c}) ' over ' ...
        num2str(nNeg+nPos) ' with fpr=' num2str(fpr) ' and fnr=' num2str(fnr)]);
    
    avgAucC = 0;
    avgAucAlt = 0;
    for i=1:nRuns
        avgAucC = avgAucC + aucInt{i}(c);
        if exist('labelsBoth','var')
            avgAucAlt = avgAucAlt + aucIntAlt{i}(c);
        end
    end
    avgAuc(c) = avgAucC ./ nRuns;
%     avgAucAlt ./ nRuns;
end
for c=1:length(combMethods)
    %calculate brier score (mean square error)
    r_anal = resultsAllFolds(c,testIdxsAllFolds);
    l_anal = labelsAllFolds(testIdxsAllFolds) - 1;

    brier(c) = sum((r_anal - l_anal).^2, 2) ./ size(l_anal,2);
end

H = zeros(length(combMethods),1);
for c=1:length(combMethods)
    r_anal = resultsAllFolds(c,testIdxsAllFolds);
    r_anal(r_anal==1) = [];
    r_anal(r_anal==0) = [];
    H(c) = -sum(r_anal.*log(r_anal) + (1-r_anal).*log(1-r_anal), 2);
end

avgAuc

