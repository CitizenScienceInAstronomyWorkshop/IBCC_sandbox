%Run this code from this directory
%clear
import settings.*

% vbIbccPaper.exp2seq;
bccSettingsDynamicsExp3;
expSettings.definiteSwitch = [42 47 52 57];

testDataFile = sprintf('%s%s_test_data.mat', expSettings.getDataDir(), ...
    expSettings.expLabel);
labelsFile = sprintf('%s%s_labels.mat', expSettings.getDataDir(), ...
    expSettings.expLabel);

if exist(testDataFile, 'file')
    labelledTestData = dlmread(testDataFile);
    labels = labelledTestData(:, expSettings.nSensors()+1);
elseif exist(labelsFile, 'file')
    labels = dlmread(labelsFile);
else
    display('cannot find labels');
end    

labels = reshape(labels, expSettings.nSamples, expSettings.nDatasets)';

combMethodsFile = sprintf('%scombMethods-%s', expSettings.getCombinerDir(), expSettings.multCombTestFile);
if exist(combMethodsFile, 'file')
    load(combMethodsFile);
else
    combMethods = {
        combiners.bcc.Ibcc.shortLabel ...
        }; 
end
% 
% combMethods = {
%     combiners.MeanDecision.shortLabel,...
%     combiners.bcc.IbccVb.shortLabel,...
%     combiners.bcc.IbccVb.shortLabelSeq,...
%     };

%loads combinedPostKnowns
load(sprintf('%s%s', expSettings.getCombinerDir(), expSettings.multCombTestFile));

agentsFile = sprintf('%s%s_multTestAgents.mat', ...
        expSettings.getDataDir(), expSettings.expLabel);
if exist(agentsFile, 'file')
    load(agentsFile);
else
    agents = cell(1, expSettings.nDatasets);
end  

basePostFile = sprintf('%s%s_base_post.mat', ...
        expSettings.getDataDir(), expSettings.expLabel);
basePost = dlmread(basePostFile);
basePost = reshape(basePost, size(basePost, 1), ...
    expSettings.nSamples, expSettings.nDatasets);    

if size(labels,1) < expSettings.nDatasets
    labels = reshape(labels, expSettings.nDatasets, expSettings.nSamples);
end

allLabels = labels;

%change the dataset here
d = 10;

datasetStart = (d-1)*expSettings.nSamples+1;
datasetEnd = d*expSettings.nSamples;
%datasetEnd = datasetStart + 199;
dataset = labelledTestData( datasetStart:datasetEnd, :);  
labels = allLabels(d, :);%1:200);

%expSettings.nSamples = 200;
%expSettings.iPropKnown = 4;

runner = ExpRunner(expSettings, [], dataset, labels+1, ...
    basePost(:, :, d), agents{d});

%change  N_knowns here.
i=1;
combinedPostRepeats = combinedPostKnowns{d, i};

runner.expSettings.iPropKnown = i;

r = 1;
drawGraph = true;
runNonDetOnly = false;

combinedPost = runner.runCombineAll(false, 'dontTouch', ...
        drawGraph, combMethods, false, runNonDetOnly);

%combinedPostR = combinedPostRepeats(r, :, :);
%combinedPost = zeros(size(combinedPostR, 3), size(combinedPostR, 2));
%for c=1:size(combinedPostR, 3)
%    combinedPost(c, :) = combinedPostR(:, :, c);
%end
% 
% chunkSize = 200;
% nPortions = expSettings.nSamples/chunkSize;
% portions = chunkSize .* (1:nPortions);

% for p=portions
%     data = combinedPost(:, p-chunkSize+1:p);
%     chunkLabels = labels(p-chunkSize+1:p);
%     
%     graph = graphs.ClassifierPerformanceGraph(expSettings, ...
%         {}, 1, combMethods, []);
%     graph.compare = false;
%     graph.drawPostGraph(data, 'Output of Combination Methods');
%     graph.drawErrorGraph(chunkLabels, data, 'Absolute Error of Combination Methods');
% end