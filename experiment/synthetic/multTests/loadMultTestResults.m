function [ labels, combinedPostKnowns, combMethods, basePost ] = loadMultTestResults( expSettings )
%LOADMULTTESTRESULTS Summary of this function goes here
%   Detailed explanation goes here
    testDataFile = sprintf('%s%s_test_data.mat', expSettings.getDataDir(), ...
        expSettings.expLabel);
    labelsFile = sprintf('%s%s_labels.mat', expSettings.getDataDir(), ...
        expSettings.expLabel);
    basePostFile = sprintf('%s%s_base_post.mat', ...
        expSettings.getDataDir(), expSettings.expLabel);
    if exist(basePostFile, 'file')
        basePost = dlmread(basePostFile);
        basePost = reshape(basePost, size(basePost, 1), ...
            expSettings.nSamples, expSettings.nDatasets);  
    else
        display('cannot find base classifier outputs');
    end
    if exist(testDataFile, 'file')
        labelledTestData = dlmread(testDataFile);
        labels = labelledTestData(:, expSettings.nSensors()+1);
    elseif exist(labelsFile, 'file')
        labels = dlmread(labelsFile);
    else
        display(sprintf('cannot find labels at %s', labelsFile));
    end    

    labels = reshape(labels(1:expSettings.nSamples*expSettings.nDatasets), expSettings.nSamples, expSettings.nDatasets)';

    combMethodsFile = sprintf('%scombMethods-%s', expSettings.getCombinerDir(), expSettings.multCombTestFile);
    if exist(combMethodsFile, 'file')
        load(combMethodsFile);
    else
        combMethods = {
            combiners.bcc.Ibcc.shortLabel ...
            }; 
    end

    load(sprintf('%s%s', expSettings.getCombinerDir(), expSettings.multCombTestFile));
end

