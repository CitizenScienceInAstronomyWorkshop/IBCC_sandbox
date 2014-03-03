%Run this code from this directory
%clear
import settings.*

%load('/homes/49/edwin/results/agent_clustering_meam-magnitude_ordered_relabelling/data/5.mat');

bccLambdaSettings;

if ~exist(expSettings.getCombinerDir(), 'file')
    mkdir(expSettings.getCombinerDir());
end

expSettings.noOverwriteData = false;
%pick the combination methods to use.
combMethods = {
    combiners.Ibcc.shortLabel ...
%    combiners.CIbcc.shortLabel...
    };

nCombiners = length(combMethods);

nTotalSamples = expSettings.nDatasets * expSettings.nSamples;

testDataFile = sprintf('%s%s_test_data.mat', ...
        expSettings.getDataDir(), expSettings.expLabel);

if exist(expSettings.getDataDir(), 'file') && exist(testDataFile, 'file')
    labelledTestData = dlmread(testDataFile);
    
    labels = labelledTestData( ...
        :, ...
        expSettings.nSensors()+1);    
else
    [labelledTrainingData labelledTestData labels] = ...
        datageneration.generateData(expSettings, nTotalSamples);
end

%uncomment this to prevent overwriting of incomplete data
%if ~exist('combinedPostLambda', 'var')
    expSettings.combinedPostLambda = cell(expSettings.nDatasets, ...
        length(expSettings.lambdaSym), length(expSettings.lambdaMag));
    combinedPostLambda = expSettings.combinedPostLambda;
%end

%runner.expSettings.lambdaMag = expSettings.lambdaMag;
%runner.expSettings.lambdaSym = expSettings.lambdaSym;
%test different values of lambda

runner.expSettings.iPropKnown = 1;

for d=1:expSettings.nDatasets
    d    
    %change datasets
    newDataset = true;
    
    dataset = labelledTestData( (d-1)*expSettings.nSamples+1: ...
        d*expSettings.nSamples, :); 
    labels = dataset(:, expSettings.nSensors()+1);
    
    runner = ExpRunner(expSettings, [], dataset, ...
        labels);
    
    for i=1:length(expSettings.lambdaSym)
        i
        
        runner.expSettings.iLambdaSym = i;
        
        for j=1:length(expSettings.lambdaMag)
            j
            
            runner.expSettings.iLambdaMag = j;
            results = zeros(expSettings.nRepeats, expSettings.nSamples);

            for r=1:expSettings.nRepeats
                r
                if i==1 && j==1 && r==1
                    results(r, :) = runner.runCombineAll(false, 'newAgents', false, combMethods);
                else
                    results(r, :) = runner.runCombineAll(false, 'dontTouch', false, combMethods);
                end
            end
        
            combinedPostLambda(d, i, j) = {results};
        end
    end
    save(sprintf('%s%s', expSettings.getCombinerDir(), expSettings.multCombTestFile), 'combinedPostLambda');

end
save(sprintf('%s%s', expSettings.getCombinerDir(), expSettings.multCombTestFile), 'combinedPostLambda');

display('done');