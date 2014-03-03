if ~exist('cyTrainAll', 'var')
    trainingFile = '/homes/49/edwin/matlab/combination/data/cyber/train.txt';
    [cyTrainAll, trainingLabelsAll] = loadCyberMatrix(trainingFile);
    testFile = '/homes/49/edwin/matlab/combination/data/cyber/test.txt';
    [cyTestAll, testLabelsAll] = loadCyberMatrix(testFile);
end
cyTrain = cyTrainAll;
cyTest = cyTestAll;
trainingLabels = trainingLabelsAll;
testLabels = testLabelsAll;

subSampleSize = 500; subSampleNAgents = 50;
cyTrain = cyTrainAll(1:subSampleNAgents, 1:subSampleSize);
trainingLabels = trainingLabelsAll(1:subSampleSize);
cyTest = cyTestAll(1:subSampleNAgents, 1:subSampleSize);
testLabels = testLabelsAll(1:subSampleSize);

% % drawPiRatings
% selectedBase = find(igRatings>2);
% cyTrain = cyTrainAll(selectedBase, :);
% cyTest = cyTestAll(selectedBase, :);
% trainingLabels = trainingLabelsAll;
% testLabels = testLabelsAll;
% subsAgents = find(sum(cyTrain, 2)~=0 & sum(cyTest,2)~=0);
% cyTrain = cyTrain(subsAgents,:);
% cyTest = cyTest(subsAgents,:);

% cyTest = sparse(size(cyTestAll,1), size(cyTestAll,2));

settings.cyber.cyber1;
expSettings.nSamples = size(cyTrain,2)+size(cyTest,2);
runner = ExpRunner(expSettings);
            
runner.noScore = -1;
runner.maxScore = 1;
runner.minScore = 0;
runner.voteThreshold = 0.5;
runner.nScores = 2;
runner.nClasses = 8;
runner.Alpha0 = [0.5 0.5; 0.5 0.5; 0.5 0.5; 0.5 0.5; 0.5 0.5; 0.5 0.5; 0.5 0.5; 0.5 0.5];

runner.scoreMap = [];
runner.noScore = -1;

drawGraphs = false;
%pick the combination methods to use.
combMethods = {...
    combiners.MeanDecision.shortLabel,...  %this will produce utter drivel for this application as features are binary, class is 8-ary and mean will always be near 0.
...%     combiners.DlrCombiner.shortLabel,...
    combiners.bcc.IbccVb.shortLabel,...   
...%     combiners.bcc.DynIbccVb.shortLabelInitAll,...
...%     combiners.bcc.IbccVb.shortLabelSeq, ...      
...%     combiners.bcc.Ibcc.shortLabel, ...
...%     combiners.weighted.WeightedSum.shortLabel,...
    };
sortResults = false;
display('starting runner...');

%display the time now
datestr(now)

%recognise that scores of zero are actually valid decisions, not no-scores
runner.combiner{2}.zeroScoreMap = 2;

testIdxs = (length(trainingLabels)+1):length(trainingLabels)+length(testLabels);

%Repeat with new feature selected until ....
allBaseData = [cyTrain cyTest];
initBaseData = [cyTrain -ones(size(cyTest,1), size(cyTest,2))];
Ass = initBaseData ~= -1; Ass = Ass';
state = {runner.Alpha0, expSettings.nu, [], Ass};

nRounds = 100;
aucs = zeros(runner.nClasses, nRounds);

for f=1:nRounds
    
    Ass = state{4};    
    currentBaseData = allBaseData; 
    vectorData = cell(1,3);
    idxs = find(Ass'~=0);
    [vectorData{1}, vectorData{2}] = ind2sub(size(currentBaseData), idxs);
    vectorData{3} = currentBaseData(idxs);
    runner.baseScoreSet = vectorData;
    display(['Number of assignments: ', num2str(sum(sum(Ass))) ]);
    
    runner.trainingLabels = [trainingLabels zeros(1, length(testLabels))]';
    runner.labels = [trainingLabels testLabels];    
    N = length(runner.labels);
    
    if f==1
        [results, agentRatings] = runner.runCombineAll(false, 'dontTouch', drawGraphs, combMethods, sortResults);
        lnPi = runner.combiner{2}.lnPi;
        lnK = runner.combiner{2}.lnK;
        ET = runner.combiner{2}.combinedPost;
    else
        [lnPi, lnK, ET, nIt] = runner.combiner{2}.iterate(1, lnPi, lnK, vectorData, size(allBaseData,2), ET);
    end
    %display the time now
    datestr(now)    
    
    testResults = cell(size(results,1),size(results,2));
    for j=1:runner.nClasses
        testResults{j} = results{j}(:, testIdxs);
    end

    auc = graphs.ClassifierPerformanceGraph.drawRoc(testResults, testLabels-1, combMethods, false, true);
    auc(:,2)
    aucs(:, f) = auc(:, 2);
%     close all;
    %%%%% PURE GREEDY (NO LOOK-AHEAD)
    K = size(cyTrain, 1);
    A = 1;  %number of allocations per user
    
    Alpha = runner.combiner{2}.Alpha;
    Nu = runner.combiner{2}.Nu;
    if runner.nClasses==2
        ET = [1-results; results];
    else
        ET = zeros(runner.nClasses, length(runner.labels));
        for j=1:runner.nClasses
            ET(j,:) = results{j}(2, :);
        end
    end
    
    state = {Alpha, Nu, ET, Ass};
    %the information gain of assigning task i to classifier k
    [chosenTask, chosenUser, IG] = greedyInfoGain(state, N, K, A, 100);
    idxs = sub2ind(size(state{4}), chosenTask, chosenUser);
    state{4}(idxs) = state{4}(idxs)+1;
end

auc = graphs.ClassifierPerformanceGraph.drawRoc(testResults, testLabels-1, combMethods, false); 
