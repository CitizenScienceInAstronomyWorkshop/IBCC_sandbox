%separate training and test phases

trainingFile = '/homes/49/edwin/matlab/combination/data/cyber/train.txt';
[cyTrain, trainingLabels] = loadCyberMatrix(trainingFile);
% 
% subSampleSize = 300;
% cyTrain = cyTrain(1:subSampleSize,1:subSampleSize);
% trainingLabels = trainingLabels(1:subSampleSize);
% cyTest = cyTest(1:subSampleSize,1:subSampleSize);
% testLabels = testLabels(1:subSampleSize);

settings.cyber.cyber1;


%TRAINING PHASE %%%%%%%%%%%%%%%%%%%%%%%%

expSettings.nSamples = size(cyTrain,2);
runner = ExpRunner(expSettings);
            
runner.noScore = -1;
runner.maxScore = 1;
runner.minScore = 0;
runner.voteThreshold = 0.5;
runner.nScores = 2;
runner.nClasses = 8;
runner.Alpha0 = [0.5 0.5; 0.5 0.5; 0.5 0.5; 0.5 0.5; 0.5 0.5; 0.5 0.5; 0.5 0.5; 0.5 0.5];

runner.basePost = cyTrain;

runner.trainingLabels = trainingLabels';
runner.labels = trainingLabels;
runner.scoreMap = [];
runner.noScore = -1;

drawGraphs = false;
%pick the combination methods to use.
combMethods = {...
    combiners.MeanDecision.shortLabel,...  %this will produce utter drivel for this application as features are binary, class is 8-ary and mean will always be near 0.
...%     combiners.DlrCombiner.shortLabel,...
    combiners.bcc.IbccVb.shortLabel,...   
...%    combiners.bcc.DynIbccVb.shortLabelInitAll,...
...%     combiners.bcc.IbccVb.shortLabelSeq, ...      
...%     combiners.bcc.Ibcc.shortLabel, ...
...%     combiners.weighted.WeightedSum.shortLabel,...
    };
sortResults = false;
display('starting runner...');

%display the time now
datestr(now)

runner.runCombineAll(false, 'dontTouch', drawGraphs, combMethods, sortResults); 

%TEST PHASE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%display the time now
datestr(now)

testFile = '/homes/49/edwin/matlab/combination/data/cyber/test.txt';
[cyTest, testLabels] = loadCyberMatrix(testFile);
% 
% subSampleSize = 300;
% cyTrain = cyTrain(1:subSampleSize,1:subSampleSize);
% trainingLabels = trainingLabels(1:subSampleSize);
% cyTest = cyTest(1:subSampleSize,1:subSampleSize);
% testLabels = testLabels(1:subSampleSize);

expSettings.nSamples = size(cyTest,2);
expSettings.nu{1} = runner.combiner{2}.Nu;
expSettings.iNu = 1;
runner2 = ExpRunner(expSettings);
runner2.Alpha0 = runner.combiner{2}.Alpha;           
runner2.noScore = -1;
runner2.maxScore = 1;
runner2.minScore = 0;
runner2.voteThreshold = 0.5;
runner2.nScores = 2;
runner2.nClasses = 8;
runner2.scoreMap = [];
runner2.noScore = -1;
runner2.trainingLabels = zeros(1, length(testLabels));
runner2.labels = testLabels;
runner2.basePost = cyTest;

[results, agentRatings] = runner2.runCombineAll(false, 'dontTouch', drawGraphs, combMethods, sortResults); 
testIdxs = 1:100;
testResults = cell(size(results,1),size(results,2));
for j=1:runner2.nClasses
    testResults{j} = results{j}(:, testIdxs);
end
%display the time now
datestr(now)

% auc = graphs.ClassifierPerformanceGraph.drawRoc(testResults, testLabels-1, combMethods, false); 
% 
% resultsDir = '/homes/49/edwin/matlab/combination/results/cyber/exp1/';
% h = get(0,'children');
% h = sort(h);
% for i=1:length(h)
%      saveas(h(i), [[resultsDir 'roc'] num2str(i)], 'png');
%      saveas(h(i), [[resultsDir 'roc'] num2str(i)], 'fig');
% end
