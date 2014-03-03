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
fpr_comb_vb = zeros(1,nRounds);

allLabels = [trainingLabels; testLabels];
currentLabels = zeros(size(allLabels)); currentLabels(1,:) = allLabels(1,:);
currentData = [cyTrain; cyTest];
unlabTrainIdxs = 2:size(trainingLabels,1);

if ~exist('nRequests','var')
    nRequests = 10; %this is currently ignored and fixed anyway
end
if ~exist('nToTry','var')
    nToTry = 50;
end
nAgents = size(cyTrain,2);
Alpha0 = ones(nClasses, 2, nAgents);
Nu0 = ones(nClasses,1);

nFaults = 0;

for f=1:nRounds
    
    display(['number of rounds: ' num2str(f)]);
           
    currentTest = [cyTest; cyTrain(unlabTrainIdxs,:)];

    nLabels = size(currentLabels, 1);

    [val,targets] = max([currentLabels ones(size(currentLabels,1),1)], [], 2);
    targets(targets>nClasses) = 0;
    
    combiner = combiners.bcc.IbccVb(nAgents, 1, targets, [], nClasses, 2, Nu0, Alpha0);
    baseData = combiner.correctBaseData(currentData');
    [combinedPost, agentRatings] = combiner.combineDecisions(baseData);
    
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
    fpr_comb_f = (log(fpr./(1-fpr)) + log(fnr./(1-fnr)))./2;
    fpr_comb_f = exp(fpr_comb_f) ./ (1+exp(fpr_comb_f))
    fpr_comb_vb(f) = fpr_comb_f;
    
%     auc = graphs.ClassifierPerformanceGraph.drawRoc(testResults, testLabVec-1, {'pseduo ibcc'}, false, true);
%     auc
%     aucs_knockon1(:, f) = auc;
%     save('/homes/49/edwin/matlab/combination/data/active_label/cyber/combinedIG.mat','aucs_combined1');

    Ptrain = P(size(cyTest,1)+1:end,:);
    
    samplesToRequest = unlabTrainIdxs(labelSelectApproxIG( nToTry, nRequests, ...
        Ptrain, currentTest(size(cyTest,1)+1:end,:), Alpha, ...
        currentLabels(size(cyTest,1)+1:end,:), nClasses ));
    
    display(num2str(samplesToRequest));
    
    currentLabels(samplesToRequest,:) = allLabels(samplesToRequest,:);
    for s=samplesToRequest
        unlabTrainIdxs(unlabTrainIdxs==s) = [];
    end    
end

auc = graphs.ClassifierPerformanceGraph.drawRoc(testResults, testLabels-1, {'pseduo ibcc'}, false); 
