%Calculate the knock-on IG of learning about one label on other labels.
%Also approximate the best set of labels to learn by greedily finding sets
%of documents that are dissimilar and also have high IG.

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

nRounds = 1000;
aucs_knockon1 = zeros(nClasses, nRounds);
lam_cov = zeros(1,nRounds);

currentLabels = trainingLabels(1,:); %randomly select a first label
currentTrain = cyTrain(1,:);
unlabTrainIdxs = 2:size(trainingLabels,1);

nRequests = 3; %this is currently ignored and fixed anyway
nToTry = 1;
nAgents = size(cyTrain,2);
Alpha0 = ones(nClasses, 2, nAgents);

nFaults = 0;

for f=1:nRounds
    
    display(['number of rounds: ' num2str(f)]);
           
    currentTest = [cyTest; cyTrain];
    currentTestLabels = [testLabels; trainingLabels];

    nLabels = size(currentLabels, 1);
    
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
    lam_cov_f = (log(fpr./(1-fpr)) + log(fnr./(1-fnr)))./2;
    lam_cov_f = exp(lam_cov_f) ./ (1+exp(lam_cov_f))
    lam_cov(f) = lam_cov_f;
    
    Ptrain = P(size(cyTest,1)+1:end,:);
    PtrUnlab = Ptrain(unlabTrainIdxs,:);

    maxCert = max(PtrUnlab, [], 2);
    [leastCert, leastCertIdx] = min(maxCert,[],1);
   
    %%%
    
    samplesToRequest = unlabTrainIdxs( stevesRoutine(X(unlabTrainIdxs,:),leastCertIdx) );

    %%%
    
    display(num2str(samplesToRequest));
    
    currentTrain = [currentTrain; cyTrain(samplesToRequest, :)];
    currentLabels = [currentLabels; trainingLabels(samplesToRequest,:)];
    for s=samplesToRequest
        unlabTrainIdxs(unlabTrainIdxs==s) = [];
    end    
end

auc = graphs.ClassifierPerformanceGraph.drawRoc(testResults, testLabels-1, {'pseduo ibcc'}, false); 
