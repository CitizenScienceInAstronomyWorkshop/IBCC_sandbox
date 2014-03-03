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
cyTrain = cyTrainAll(5001:subSampleSize+5001, 1:subSampleNAgents);
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
fpr_comb2 = zeros(1,nRounds);

currentLabels = trainingLabels(1,:); %randomly select a first label
currentTrain = cyTrain(1,:);
unlabTrainIdxs = 2:size(trainingLabels,1);

if ~exist('nRequests','var')
    nRequests = 10; %this is currently ignored and fixed anyway
end
if ~exist('nToTry','var')
    nToTry = 50;
end
nAgents = size(cyTrain,2);
Alpha0 = ones(nClasses, 2, nAgents);

nFaults = 0;

for f=1:nRounds
    
    display(['number of rounds: ' num2str(f)]);
           
    currentTest = [cyTest; cyTrain(unlabTrainIdxs,:)];
    currentTestLabels = [testLabels; trainingLabels(unlabTrainIdxs,:)];

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
    fpr_comb_f = (log(fpr./(1-fpr)) + log(fnr./(1-fnr)))./2;
    fpr_comb_f = exp(fpr_comb_f) ./ (1+exp(fpr_comb_f))
    fpr_comb2(f) = fpr_comb_f;
    
%     auc = graphs.ClassifierPerformanceGraph.drawRoc(testResults, testLabVec-1, {'pseduo ibcc'}, false, true);
%     auc
%     aucs_knockon1(:, f) = auc;
%     save('/homes/49/edwin/matlab/combination/data/active_label/cyber/combinedIG.mat','aucs_combined1');

    %%%%% PURE GREEDY (NO LOOK-AHEAD)
    %select a label from the **training data** that we haven't yet seen
    Ptrain = P(size(cyTest,1)+1:end,:);

    maxCert = max(Ptrain, [], 2);
    [leastCert, leastCertIdx] = sort(maxCert,1,'ascend');
    
    samplesToTry = unlabTrainIdxs(leastCertIdx(1:nToTry));
    PtoTry = Ptrain(leastCertIdx(1:nToTry), :);
    
    
    IG = zeros(size(P,1),length(samplesToTry));
    Pafter = zeros(size(P,1), size(P,2), nClasses);
    
    nScores = size(Alpha,2);
    
    for i=1:length(samplesToTry)
        s = samplesToTry(i); %index to training samples
        P_s = Ptrain(leastCertIdx(i),:); %probability of the training sample class according to last run of classifier
                
        otherUnlabTrainIdxs = unlabTrainIdxs(unlabTrainIdxs~=s);
        sampleTestData = [cyTest; cyTrain(otherUnlabTrainIdxs,:)];
               
        [AlphaMatList, lnAlphaMatList, AlphaMatTotalsList] = getAlphaInCmatForm(Alpha,currentTest+1);

        for j=1:nClasses
            Nu = sum(currentLabels,1)+1; Nu(j) = Nu(j) + 1;
            NuSum = sum(Nu);
            Kappa = Nu/NuSum;
            [IGj, eBefore, eAfter, Pafterj] = objectInfoGain(Kappa, P, P_s, currentTest+1, Alpha, ...
                j, cyTrain(s,:)+1, AlphaMatList, lnAlphaMatList, AlphaMatTotalsList);
            IG(:,i) = IG(:,i) + P_s(j).*IGj;
            Pafter(:,:,j) = Pafterj;
        end
        ownIdx = size(cyTest,1) + leastCertIdx(i);
        IG(ownIdx,i) = IG(ownIdx,i) - sum(sum(Pafter(ownIdx, :, :),3).*log(P_s));        
%         if ~isempty(find(IG(:,i)<0,1))
%             nFaults = nFaults + 1;
%             save('/homes/49/edwin/matlab/combination/data/active_label/cyber/nFaults.mat','nFaults');
%         end        
    end
    
    totalIG = sum(IG,1);
    if ~isempty(find(totalIG<-0.1,1))
        display('Error: negative IG!');
    end
    testSampleIdxs = size(cyTest,1) + leastCertIdx(1:nToTry);
    totalValue = zeros(1,nToTry);
    sampleSet = zeros(3,nToTry);
    for i=1:length(samplesToTry)
        
        selected = i;
        
        for r=1:nRequests-1
            jointIGEst = sum(IG(testSampleIdxs,selected),2) + max(IT(testSampleIdxs,selected),[],2);
            [MInext MInextIdx] = min(jointIGEst);
            selected = [selected; MInextIdx];
        end
        
        sampleSet(:,i) = selected;
        totalValue(i) = MInext;
        
%         [MIj minMIIdx] = min(IG(testSampleIdxs,i));
%         j = minMIIdx;
%         jointMI = IG(testSampleIdxs,i)+IG(testSampleIdxs,j) +abs(IG(testSampleIdxs,i)-IG(testSampleIdxs,j));%+ max(IG(testSampleIdxs,i)+IG(testSampleIdxs,j));
%         [MIk minMIIdx] = min(jointMI);
%         k = minMIIdx;
%         totalValue(i) = sum(IG(:,i)+IG(:,j)+IG(:,k)) + sum(max([IG(:,i) IG(:,j) IG(:,k)], [], 2));%totalIG(i) + totalIG(j)/(1+MIj) + totalIG(k)/(1+MIk);
%         
%         sampleSet(1,i) = i;
%         sampleSet(2,i) = j;
%         sampleSet(3,i) = k;
        
        display([num2str(i) ', ' num2str(sampleSet(:,i)') ', ' num2str(sum(totalIG(i))) ', ' num2str(sum(totalValue(i))) ', ' num2str(leastCert(i))]);        
    end    
    
    [maxVal setToRequest] = max(totalValue);
    sampleToRequest = samplesToTry(sampleSet(:,setToRequest));
    display(num2str(sampleToRequest));
    
    currentTrain = [currentTrain; cyTrain(sampleToRequest, :)];
    currentLabels = [currentLabels; trainingLabels(sampleToRequest,:)];
    for s=sampleToRequest
        unlabTrainIdxs(unlabTrainIdxs==s) = [];
    end    
end

auc = graphs.ClassifierPerformanceGraph.drawRoc(testResults, testLabels-1, {'pseduo ibcc'}, false); 
