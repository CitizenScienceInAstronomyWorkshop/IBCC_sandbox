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

nRounds = 60;
aucs = zeros(nClasses, nRounds);
lam_uncert = zeros(1,nRounds);
currentLabels = trainingLabels(1,:); %randomly select a first label
currentTrain = cyTrain(1,:);
trainIdxs = 1:size(trainingLabels,1);

nRequests = 3;
nRepeats = 1; %repeatedly request the same label because if might be faulty - get a majority opinion
chosenOnes = [];

for f=1:nRounds
       
    currentTest = [cyTest; cyTrain];
    currentTestLabels = [testLabels; trainingLabels];
    
    [P_class,P_feat_class,Alpha,Nu] = multi_classifier(currentTrain,currentLabels);
    P = multi_classify(P_class,P_feat_class,currentTest);

    %display the time now
    datestr(now)    
    
    testResults = cell(size(P,2),1);
    for j=1:nClasses
        testResults{j} = P(:, j)';
    end

%     auc = graphs.ClassifierPerformanceGraph.drawRoc(testResults, testLabVec-1, {'pseduo ibcc'}, false, true);
%     auc
%     aucs(:, f) = auc;
    [maxCert maxClass] = max(P,[],2);
    binaryRes = zeros(size(P));
    binaryRes(sub2ind(size(P),1:size(P,1),maxClass')) = 1;
    fpr = sum(sum((binaryRes-currentTestLabels)>0)) ./ sum(sum(currentTestLabels==0));
    fnr = sum(sum((binaryRes-currentTestLabels)<0)) ./ sum(sum(currentTestLabels==1));
    lam_uncert_f = (log(fpr./(1-fpr)) + log(fnr./(1-fnr)))./2;
    lam_uncert_f = exp(lam_uncert_f) ./ (1+exp(lam_uncert_f))
    
    for rep=0:(nRepeats-1)
        lam_uncert(f*nRepeats - rep) = lam_uncert_f;
    end

    %%%%% PURE GREEDY (NO LOOK-AHEAD)
    %select a label from the **training data** that we haven't yet seen
    for r=1:nRequests
        Ptr = P(size(cyTest,1)+1:end,:);
    
        notChosenOnes = trainIdxs(~ismember(trainIdxs,chosenOnes));
        Ptr = Ptr(notChosenOnes,:);
        
        
%         while ismember(leastCertIdx, chosenOnes)
%             leastCertIdx = randi(length(trainIdxs),1);
%         end        
        
        maxCert = max(Ptr, [], 2);
        [leastCert, leastCertIdx] = min(maxCert);
        
        leastCertIdx = notChosenOnes(leastCertIdx);
    
        for rep=1:nRepeats
        currentTrain = [currentTrain; cyTrain(trainIdxs(leastCertIdx), :)];
        
        newLabel = trainingLabels(trainIdxs(leastCertIdx),:);
%         pMutate = 0.3;
%         mutate = rand(1,1);
%         if mutate < pMutate
%             c = find(newLabel==1);
%             newLabel(c) = 0;
%             if c<nClasses
%                 newLabel(c+1) = 1;
%             else
%                 newLabel(1) = 1;
%             end
%         end
        
        currentLabels = [currentLabels; newLabel ];
        chosenOnes = [chosenOnes trainIdxs(leastCertIdx)];
        end
    end
    
end

% auc = graphs.ClassifierPerformanceGraph.drawRoc(testResults, testLabels-1, {'pseduo ibcc'}, false); 
