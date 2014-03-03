if ~exist('cyTrainAll', 'var')
    trainingFile = '/homes/49/edwin/matlab/data/cyber/train.txt';
    [cyTrainAll, trainingLabelsAll] = loadCyberMatrix(trainingFile);
    testFile = '/homes/49/edwin/matlab/data/cyber/test.txt';
    [cyTestAll, testLabelsAll] = loadCyberMatrix(testFile);
end
% cyTrain = cyTrainAll;
% cyTest = cyTestAll;
% trainingLabels = trainingLabelsAll;
% testLabels = testLabelsAll;

% subSampleSize = 1000;
% cyTrain = cyTrainAll(1:subSampleSize,1:subSampleSize);
% trainingLabels = trainingLabelsAll(1:subSampleSize);
% cyTest = cyTestAll(1:subSampleSize,1:subSampleSize);
% testLabels = testLabelsAll(1:subSampleSize);

% drawPiRatings
if exist('igRatings','var')
    selectedBase = find(igRatings>0.0020);
    cyTrain = cyTrainAll(selectedBase, :);
    cyTest = cyTestAll(selectedBase, :);
    display('Selecting only a subset of base classifiers with IG ratings > 0.0020');
else
    cyTrain = cyTrainAll;
    cyTest = cyTestAll;
end
trainingLabels = trainingLabelsAll;
testLabels = testLabelsAll;

subsAgents = find(sum(cyTrain, 2)~=0 & sum(cyTest,2)~=0);
cyTrain = cyTrain(subsAgents,:);
cyTest = cyTest(subsAgents,:);%zeros(size(cyTrain,1),0);%

settings.cyber.cyber1;
expSet.nSamples = size(cyTrain,2)+size(cyTest,2);

testIdxs = (length(trainingLabels)+1):length(trainingLabels)+length(testLabels);
runnerBasePost = [cyTrain cyTest];

runnerTrLabels = [trainingLabels zeros(1, length(testLabels))]';
labels = [trainingLabels testLabels] + 1;

runner = ExpRunner(expSet, [], bccSet, [], [], labels, runnerBasePost, [], runnerTrLabels);

drawGraphs = false;
%pick the combination methods to use.
combMethods = {...
...%     combiners.MeanDecision.shortLabel,...  %this will produce utter drivel for this application as features are binary, class is 8-ary and mean will always be near 0.
         combiners.LrCombiner.shortLabel,...
...%     combiners.DlrCombiner.shortLabel,...
         combiners.bcc.IbccVb.shortLabel,...   
...%     combiners.bcc.IbccSampling.shortLabel, ...
    };
sortResults = false;
display('starting runner...');

%display the time now
datestr(now)

%recognise that scores of zero are actually valid decisions, not no-scores
runner.combiner{2}.zeroScoreMap = 2;

[results, agentRatings] = runner.runCombineAll(false, 'dontTouch', drawGraphs, combMethods, sortResults); 
testResults = cell(size(results,1),size(results,2));
for j=1:expSet.nClasses
    testResults{j} = results{j}(:, testIdxs);
end
%display the time now
datestr(now)

auc = graphs.ClassifierPerformanceGraph.drawRoc(testResults, testLabels-1, combMethods, false); 

resultsDir = '/homes/49/edwin/matlab/results/cyber/exp1_topfeaturesonly/';
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end
h = get(0,'children');
h = sort(h);
for i=1:length(h)
     saveas(h(i), [[resultsDir 'roc'] num2str(i)], 'png');
     saveas(h(i), [[resultsDir 'roc'] num2str(i)], 'fig');
end

% save('/homes/49/edwin/matlab/combination/data/cyber/
