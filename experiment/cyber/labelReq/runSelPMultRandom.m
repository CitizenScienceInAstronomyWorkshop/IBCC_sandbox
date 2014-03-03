%using pseduo-count update method with a separate training phase instead of
%VB. VB could be used later once labels have been selected.

load('Cyber Data/train.mat');
trainingLabelsAll = T;
cyTrainAll = X;

load('Cyber Data/test.mat');
testLabelsAll = T;
cyTestAll = X;

nClasses = size(T,2);
nFeatures = size(X,2);

subSampleSize = 5000; subSampleNAgents = 500;
cyTrain = cyTrainAll(1:subSampleSize, 1:subSampleNAgents);
trainingLabels = trainingLabelsAll(1:subSampleSize,:);
[tmp, trainLabVec] = max(trainingLabels,[],2);

% cyTest = cyTestAll(1:subSampleSize, 1:subSampleNAgents);
% testLabels = testLabelsAll(1:subSampleSize,:);
% [tmp, testLabVec] = max(testLabels,[],2);
cyTest = zeros(0,subSampleNAgents);
testLabels = zeros(0,nClasses);

settings.cyber.cyber1;
expSettings.nSamples = size(cyTrain,2)+size(cyTest,2);

%display the time now
datestr(now)

%Repeat with new feature selected until ....
allBaseData = [cyTrain; cyTest];
initBaseData = [cyTrain; -ones(size(cyTest,1), size(cyTest,2))];
Ass = initBaseData ~= -1;
state = {ones(nClasses,2,nFeatures), expSettings.nu, [], Ass};

nRounds = 300;
aucs = zeros(nClasses, nRounds);
lam_rand3 = zeros(1,nRounds);
currentLabels = trainingLabels(1,:); %randomly select a first label
currentTrain = cyTrain(1,:);
trainIdxs = 1:size(trainingLabels,1);

nRequests = 3;

chosenOnes = [];

for f=1:nRounds
    notChosenOnes = trainIdxs(~ismember(trainIdxs,chosenOnes));
    currentTest = [cyTest; cyTrain(notChosenOnes,:)];
    currentTestLabels = [testLabels; trainingLabels(notChosenOnes,:)];
    
    [P_class,P_feat_class,Alpha,Nu] = multi_classifier(currentTrain,currentLabels);
    P = multi_classify(P_class,P_feat_class,currentTest);

    %display the time now
    datestr(now)    
    
    testResults = cell(size(P,2),1);
    for j=1:nClasses
        testResults{j} = P(:, j)';
    end
    [maxCert maxClass] = max(P,[],2);
    binaryRes = zeros(size(P));
    binaryRes(sub2ind(size(P),1:size(P,1),maxClass')) = 1;
    fpr = sum(sum((binaryRes-currentTestLabels)>0)) ./ sum(sum(currentTestLabels==0));
    fnr = sum(sum((binaryRes-currentTestLabels)<0)) ./ sum(sum(currentTestLabels==1));
    lam_rand_f = (log(fpr./(1-fpr)) + log(fnr./(1-fnr)))./2;
    lam_rand_f = exp(lam_rand_f) ./ (1+exp(lam_rand_f))
    lam_rand3(f) = lam_rand_f;
%     auc = graphs.ClassifierPerformanceGraph.drawRoc(testResults, testLabVec-1, {'pseduo ibcc'}, false, true);
%     auc
%     aucs(:, f) = auc;

    %%%%% PURE GREEDY (NO LOOK-AHEAD)
    %select a label from the **training data** that we haven't yet seen
    for r=1:nRequests
        leastCertIdx = randi(length(trainIdxs),1);
%         while ismember(leastCertIdx, chosenOnes)
%             leastCertIdx = randi(length(trainIdxs),1);
%         end
    
        currentTrain = [currentTrain; cyTrain(trainIdxs(leastCertIdx), :)];
        currentLabels = [currentLabels; trainingLabels(trainIdxs(leastCertIdx),:)];
%         chosenOnes = [chosenOnes leastCertIdx];
    end
    
end

% auc = graphs.ClassifierPerformanceGraph.drawRoc(testResults, testLabels-1, {'pseduo ibcc'}, false); 
