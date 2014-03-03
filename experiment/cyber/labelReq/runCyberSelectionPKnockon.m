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
aucs_knockon1 = zeros(nClasses, nRounds);
fpr_knockon = zeros(1,nRounds);
currentLabels = trainingLabels(1,:); %randomly select a first label
currentTrain = cyTrain(1,:);
unlabTrainIdxs = 2:size(trainingLabels,1);

nRequests = 3;
nToTry = 10;
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
    fpr_knockon_f = (log(fpr./(1-fpr)) + log(fnr./(1-fnr)))./2;
    fpr_knockon_f = exp(fpr_knockon_f) ./ (1+exp(fpr_knockon_f))
    fpr_knockon(f) = fpr_knockon_f;
%     
%     auc = graphs.ClassifierPerformanceGraph.drawRoc(testResults, testLabVec-1, {'pseduo ibcc'}, false, true);
%     auc
%     aucs_knockon1(:, f) = auc;
%     save('/homes/49/edwin/matlab/combination/data/active_label/cyber/knockonIG_3req.mat','aucs_knockon1');

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
               
%         newCTest = cell(1,3);
%         newCTest{1} = reshape(repmat((1:size(sampleTestData,2))', size(sampleTestData,1), 1), 1, numel(sampleTestData));
%         newCTest{2} = repmat(1:size(sampleTestData,1), 1, size(sampleTestData,2));
%         newCTest{3} = reshape(sampleTestData, 1, numel(sampleTestData)) + 1;    

        [AlphaMatList, lnAlphaMatList, AlphaMatTotalsList] = getAlphaInCmatForm(Alpha,currentTest+1);
%         PiIndx = sub2ind([nScores nAgents], sampleTestData, repmat(1:nAgents,size(sampleTestData,1),1));   


        for j=1:nClasses
            %newET,CTrain
%             sampleLabel = zeros(1,nClasses); sampleLabel(j) = 1;
            
            %need to adjust P_class given the observed label?
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
        if ~isempty(find(IG(:,i)<0,1))
            nFaults = nFaults + 1;
            save('/homes/49/edwin/matlab/combination/data/active_label/cyber/nFaults.mat','nFaults');
        end        
        display([num2str(i) ', ' num2str(s) ', ' num2str(sum(IG(:,i))) ', ' num2str(leastCert(i))]);
    end
    
    totalIG = sum(IG,1);
    if ~isempty(find(totalIG<-0.1,1))
        display('Error: negative IG!');
    end
    
    [sortedIG, sortedIGIdx] = sort(totalIG,2,'descend');
    sampleToRequest = samplesToTry(sortedIGIdx(1:nRequests));
    display(num2str(sampleToRequest));
    
    currentTrain = [currentTrain; cyTrain(sampleToRequest, :)];
    currentLabels = [currentLabels; trainingLabels(sampleToRequest,:)];
    for s=sampleToRequest
        unlabTrainIdxs(unlabTrainIdxs==s) = [];
    end    
end

% auc = graphs.ClassifierPerformanceGraph.drawRoc(testResults, testLabels-1, {'pseduo ibcc'}, false); 
