% !!! The base test and training data needs to be recreated for each setting!!!

%To use this script, first load the settings object and the
%combMethods object.

if ~exist('drawGraph','var')
    drawGraph = false;
end

combMethods = reshape(combMethods, 1, numel(combMethods));

save(sprintf('%scombMethods-%s', expSet.getCombinerDir(), expSet.multCombTestFile), ...
    'combMethods');

nCombiners = numel(combMethods);

%set newdata argument to false to reuse the previous data; newagents
%argument to false to reuse existing agent-sensor assignments and training.

%test different numbers of known labels

nTotalSamples = synthSet.nDatasets * expSet.nSamples;

testDataFile = sprintf('%s%s_test_data.mat', ...
        expSet.getDataDir(), expSet.expLabel);

%store the outputs of the base classifiers. Managing this at this level
%because experiment runner doesn't handle multiple datasets.
basePostFile = sprintf('%s%s_base_post.csv', ...
        expSet.getDataDir(), expSet.expLabel);
    
labelsFile = sprintf('%s%s_labels.csv', ...
        expSet.getDataDir(), expSet.expLabel);
    
agentsFile = sprintf('%s%s_multTestAgents.mat', ...
        expSet.getDataDir(), expSet.expLabel);

    
if ~exist('dataLengths', 'var') && ~exist('nSettings','var')
    nSettings = length(expSet.propKnown);
elseif ~exist('nSettings','var')
    nSettings = length(dataLengths);
    expSet.iPropKnown = 1;
end

if ~exist('labelsKnowns','var') || size(labelsKnowns,1) ~=synthSet.nDatasets || size(labelsKnowns,2) ~= nSettings
  labelsKnowns = cell(synthSet.nDatasets, nSettings);
end

if ~exist('combinedPostKnowns','var') || size(combinedPostKnowns,1) ~= synthSet.nDatasets || size(combinedPostKnowns,2) ~= nSettings
    combinedPostKnowns = cell(synthSet.nDatasets, nSettings);
    agentPostKnowns = cell(synthSet.nDatasets, nSettings);

    settingStart = 1;
    runStart = 1;
else
    settingStart = 1;
    while settingStart<size(combinedPostKnowns,2) && ~isempty(combinedPostKnowns{1,settingStart+1})
        settingStart = settingStart + 1;
    end
    runStart = 1;
end

progressMonitor = scoring.CombinerProgressMonitor([], length(combMethods));

if ~exist('rerunAgentsEachIt','var')
    rerunAgentsEachIt = false;  
end

errorSum = zeros(synthSet.nAgents, nSettings);
fnr_agents_sum = zeros(synthSet.nAgents, nSettings);
fpr_agents_sum = zeros(synthSet.nAgents, nSettings);
combinedError = zeros(length(combMethods), nSettings);
combinedErrorAll = zeros(length(combMethods), nSettings, synthSet.nDatasets);
entropy = zeros(length(combMethods), nSettings);
fpr = zeros(length(combMethods), nSettings);
fnr = zeros(length(combMethods), nSettings);
   
for i=settingStart:nSettings

    if exist('varySettingFunc','var')
        [expSet synthSet bccSet, rerunAgentsEachIt] = varySettingFunc(i);
        if ~exist('rootSuffix','var')
            rootSuffix = [];
        end     
        replacementDir = [expSet.rootDir rootSuffix 'var_' num2str(i)];
        expSet.setDirNames(replacementDir);
    elseif ~exist('dataLengths', 'var')
        expSet.iPropKnown = i;
    else
        expSet.nSamples = dataLengths(i);            
    end
    
    nTrSamplesAllAgents = synthSet.nTrainingSamples * synthSet.nAgents;
    nTotalTrSamples = synthSet.nDatasets * nTrSamplesAllAgents;    

    if ~rerunAgentsEachIt && i==1 && exist(expSet.getDataDir(), 'file') ...
            && exist(testDataFile, 'file')
        labelledTestData = dlmread(testDataFile);

        if exist(basePostFile, 'file')
            basePost = dlmread(basePostFile);
            basePost = reshape(basePost, size(basePost, 1), ...
                expSet.nSamples, synthSet.nDatasets);        
            rerunBaseAgents = false;  
        else
            basePost = zeros(expSet.nAgents, expSet.nSamples, ...
            synthSet.nDatasets);
            rerunBaseAgents = true;
        end

        if exist(agentsFile, 'file')
            load(agentsFile);
        else
            agents = cell(1, synthSet.nDatasets);
        end
    elseif rerunAgentsEachIt || i==1
        
        synthSet.nTrainingSamples = nTotalTrSamples;
        [labelledTrainingData labelledTestData] = ...
            datageneration.generateData(synthSet, expSet.expLabel, nTotalSamples, synthSet.nAgents);
        basePost = zeros(synthSet.nAgents, expSet.nSamples, synthSet.nDatasets);
        rerunBaseAgents = true;

        synthSet.nTrainingSamples = nTrSamplesAllAgents./synthSet.nAgents;
        
        agents = cell(1, synthSet.nDatasets);
    end

    runStart = 1;
    while runStart<size(combinedPostKnowns,1) && ~isempty(combinedPostKnowns{runStart,i})
        runStart = runStart + 1;
    end

    if runStart==size(combinedPostKnowns,1) && ~isempty(combinedPostKnowns{runStart,i})
        continue;
    end

    display(['Settings starting from ' num2str(settingStart) ', datasets from ' num2str(runStart)]);

    for d=runStart:synthSet.nDatasets
        %change datasets
        %newDataset = true;
        display(sprintf('d=%d', d));
        datasetStart = (d-1)*expSet.nSamples+1;
        datasetEnd = d*expSet.nSamples;
        dataset = labelledTestData( datasetStart:datasetEnd, :);  
        
        trDatasetStart = (d-1)*nTrSamplesAllAgents+1;
        trDatasetEnd = d*nTrSamplesAllAgents;
        agentTrainingSet = labelledTrainingData( trDatasetStart:trDatasetEnd, :);
        
        labels = dataset(:, synthSet.nSensors()+1)' + 1;

                
        if rerunBaseAgents || rerunAgentsEachIt
            runner = ExpRunner(expSet, synthSet, bccSet, agentTrainingSet, dataset, labels);
        else
            runner = ExpRunner(expSet, synthSet, bccSet, agentTrainingSet, dataset, labels, ...
                basePost(:, :, d), agents{d});
        end

        progressMonitor.updateLabels(labels);
        runner.progressMonitor = [];%progressMonitor;
            
%          display('what was the dataset label for?');
        %         runner.datasetLabel = sprintf('d%i_knowns%i', d, i);
        
        display(sprintf('i=%d, d=%d', i, d));   
        for r=1:expSet.nRepeats
            display(sprintf('r=%d', r));
                        
            if (rerunBaseAgents && r==1 && i==1) || rerunAgentsEachIt
                
                [~, agentRatings, ~, basePost_r, runResults, l, baseKnowns] = runner.runCombineAll(...
                    false, 'newAgents', drawGraph, combMethods);
                
                results = zeros(expSet.nRepeats, length(l), nCombiners);
                
                basePost(:, 1:size(basePost_r,2), d) = basePost_r(1:size(basePost,1), :, :);
                
                agents{d} = runner.agents';
%                 save(agentsFile, 'agents');

                basePostFile = sprintf('%s../../data/base_var%d_d%d.csv', ...
                            expSet.getDataDir(), i, d);
    
                labelsFile = sprintf('%s../../data/labels_var%d_d%d.csv', ...
                            expSet.getDataDir(), i, d);
                dlmwrite(basePostFile, basePost_r');
                dlmwrite(labelsFile, labels');
                
                for c=1:nCombiners
                    results(r, :, c) = runResults(c, :);
                end
                
                [errors fnr_agents fpr_agents] = runner.agentErrorRates();
                errorSum(:,i) = errorSum(:,i) + errors(1:size(errorSum,1));
                fnr_agents_sum(:,i) = fnr_agents_sum(:,i) + fnr_agents(1:size(fnr_agents_sum,1));
                fpr_agents_sum(:,i) = fpr_agents_sum(:,i) + fpr_agents(1:size(fpr_agents_sum,1)); 
            else       
                if r > 1
                    runNonDetOnly = true; %only re-run methods that might change
                else
                    runNonDetOnly = false;
                end
                [~, agentRatings, ~, ~, runResults, l, baseKnowns] = runner.runCombineAll(false, 'dontTouch', ...
                    drawGraph, combMethods, false, runNonDetOnly);
                
                if r==1
                    results = zeros(expSet.nRepeats, length(l), nCombiners);                
                end
                
                for c=1:nCombiners
                    %if not re-run the results will be set to -1
                    if sum(runResults(c, :)==-1) == size(runResults, 2)
                        runResults(c, :) = results(r-1, :, c);
                    end
                    
                    results(r, :, c) = runResults(c, :);
                end
            end
        end
        l_anal = repmat(l, size(runResults,1), 1);
        combinedErrorAll(:,i,d) = sum((runResults - l_anal + 1).^2, 2)./size(runResults,2);
        combinedError(:, i) = combinedError(:, i) + combinedErrorAll(:,i,d);
        
        logRR = log(runResults); logRR(isinf(logRR)) = 0;
        log1RR = log(1-runResults); log1RR(isinf(log1RR)) = 0;
        entropy(:, i) = entropy(:, i) - sum(runResults.*logRR + (1-runResults).*log1RR,2);
        posIdxs = l_anal(1,:)==2;
        nPos = sum(posIdxs,2);
        fnr(:,i) = fnr(:,i) + sum(1-round(runResults(:,posIdxs)),2) ./ nPos;
        
        negIdxs = l_anal(1,:)==1;
        nNeg = sum(negIdxs,2);
        fpr(:,i) = fpr(:,i) + sum(round(runResults(:,negIdxs)),2) ./ nNeg;
        
        labelsKnowns(d,i) = {l};
        
        agentPostKnowns(d, i) = {baseKnowns};
        combinedPostKnowns(d, i) = {results};
        
        save(sprintf('%s%s', expSet.getCombinerDir(), expSet.multCombTestFile), ...
            'combinedPostKnowns');
        combs = runner.combiner;
        save([outputDir '/experiment_' num2str(expNo) '_d' num2str(d) '_i' num2str(i) '_combiners.mat'], 'combs');        
    end  
    rerunBaseAgents = false;
end
combinedError = combinedError ./ synthSet.nDatasets;
entropy = entropy ./ synthSet.nDatasets;
fnr = fnr ./ synthSet.nDatasets;
fpr = fpr ./ synthSet.nDatasets;

errorSum = errorSum ./ synthSet.nDatasets;
bestAgentError = min(errorSum,1);%sum(errorSum,1) ./ synthSet.nAgents;
meanAgentError = sum(errorSum,1)./size(errorSum,1);

fpr_agents_sum = fpr_agents_sum ./ synthSet.nDatasets;
bestAgentFpr = min(fpr_agents_sum,1);
fnr_agents_sum = fnr_agents_sum ./ synthSet.nDatasets;
bestAgentFnr = min(fnr_agents_sum,1);

save(sprintf('%s%s', expSet.getCombinerDir(), expSet.multCombTestFile), ...
    'combinedPostKnowns', 'labelsKnowns');

save(sprintf('%s%s', expSet.getCombinerDir(), 'progressMonitor'), 'progressMonitor');

display('done');