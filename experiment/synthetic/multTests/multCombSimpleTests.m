%To use this script, first load the settings object and the
%combMethods object.

combMethods = reshape(combMethods, 1, numel(combMethods));

save(sprintf('%scombMethods-%s', expSettings.getCombinerDir(), expSettings.multCombTestFile), ...
    'combMethods');

nCombiners = numel(combMethods);

%set newdata argument to false to reuse the previous data; newagents
%argument to false to reuse existing agent-sensor assignments and training.

%test different numbers of known labels

nTotalSamples = expSettings.nDatasets * expSettings.nSamples;

%test data 
testDataFile = sprintf('%s%s_test_data.mat', ...
        expSettings.getDataDir(), expSettings.expLabel);
%just contains labels
labelsFile = sprintf('%s%s_labels.mat', ...
        expSettings.getDataDir(), expSettings.expLabel);

%fake agents: this file contains probability of correct answer for each
%fake agent.
agentsFile = sprintf('%s%s_multTestAgents.mat', ...
        expSettings.getDataDir(), expSettings.expLabel);

%the outputs of the fake agents
basePostFile = sprintf('%s%s_base_post.mat', ...
        expSettings.getDataDir(), expSettings.expLabel);
        
if exist(expSettings.getDataDir(), 'file')
    if exist(testDataFile, 'file')
        labelledTestData = dlmread(testDataFile);
        labels = labelledTestData( ...
            :, ...
            expSettings.nSensors()+1);
    elseif exist(labelsFile, 'file')
        labels = dlmread(labelsFile);
    end
    
    if exist(basePostFile, 'file')
        basePost = dlmread(basePostFile);
        basePost = reshape(basePost, size(basePost, 1), ...
            expSettings.nSamples, expSettings.nDatasets);        
        rerunBaseAgents = false;
    else
        basePost = zeros(expSettings.nAgents, expSettings.nSamples, ...
        expSettings.nDatasets);
        rerunBaseAgents = true;  
    end    
    
    if exist(agentsFile, 'file')
        %load them in case we need to re-calculate basePost
        load(agentsFile);
    else
        %experiment runner will create them, we will save them here
        agents = cell(1, expSettings.nDatasets);
    end
else
    %generate two-class labels
    labels = datageneration.generateLabels(2, nTotalSamples, expSettings);
    
    %experiment runner will create them, we will save them here
    agents = cell(1, expSettings.nDatasets);
    
    %create and run the fake agents - the results go here in basePost
    basePost = zeros(expSettings.nAgents, expSettings.nSamples, ...
        expSettings.nDatasets);
    rerunBaseAgents = true;
end

combinedPostKnowns = expSettings.combinedPostKnowns;

if ~exist('dataLengths', 'var')
    nSettings = length(expSettings.propKnown);
else
    nSettings = length(dataLengths);
    expSettings.iPropKnown = 1;
end

% labels = labels + 1;

%run combination for each dataset
for d=1:expSettings.nDatasets
    %change datasets
    %newDataset = true;
    display(sprintf('d=%d', d));
    datasetStart = (d-1)*expSettings.nSamples+1;
    datasetEnd = d*expSettings.nSamples;
      
    datasetLabels = labels(datasetStart:datasetEnd);
    
    if rerunBaseAgents
        runner = ExpRunner(expSettings, [], [], datasetLabels);
    else
        runner = ExpRunner(expSettings, [], [], datasetLabels, ...
            basePost(:, :, d));
    end
    
    if exist('Alpha0', 'var')
        runner.Alpha0 = Alpha0;
    end
        
    %run combination for each setting
    for i=1:nSettings
        %vary amount of test data or proportion of known labels?
        if ~exist('dataLengths', 'var')
            runner.expSettings.iPropKnown = i;
        else
            runner.expSettings.nSamples = dataLengths(i);
        end
        
        results = zeros(expSettings.nRepeats, expSettings.nSamples, ...
            nCombiners);
        
        display(sprintf('i=%d, d=%d', i, d));
        %run several repeats for combination methods that use stochastic methods 
        for r=1:expSettings.nRepeats
            display(sprintf('r=%d', r));
            
            if rerunBaseAgents && r==1 && i==1
                [runResults, agentRatings, l, basePost(:, :, d)] = runner.runCombineAll(...
                    false, 'newFakeAgents', false, combMethods);
                
                agents{d} = runner.agents';
                save(agentsFile, 'agents');
                dlmwrite(basePostFile, basePost);
                
                for c=1:nCombiners
                    results(r, :, c) = runResults(c, :);
                end
            else
                runResults = runner.runCombineAll(false, 'dontTouch', ...
                    false, combMethods);
                for c=1:nCombiners
                    results(r, :, c) = runResults(c, :);
                end
            end
        end
        
        combinedPostKnowns(d, i) = {results};
        
        save(sprintf('%s%s', expSettings.getCombinerDir(), expSettings.multCombTestFile), ...
            'combinedPostKnowns');
    end
end
save(sprintf('%s%s', expSettings.getCombinerDir(), expSettings.multCombTestFile), ...
    'combinedPostKnowns');

display('done');