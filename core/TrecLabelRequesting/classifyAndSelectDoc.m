function [samplesToRequest, P, binaryRes, pairsRes, origTsorted] = classifyAndSelectDoc( nToTry, nRequests, ...
    featureFile, labelFile, outputRequestFile, outputResultFile, ...
    topicDocNoPairsFile, fileMapFile, Alpha0, Nu0, ...
    writeOutput, selectMethod, classifierType, confWeight, crowdTrust, runName )

% display('classifying then selecting next document to label...');

if nargin < 15
    runName = 'Orchid2012';
end

if nargin < 11
    writeOutput = true;
end

if nargin < 12
    selectMethod = 'uncertClust';
end

if nargin < 13
    classifierType = 'two-stage';
end

if nargin < 10 || isempty(Nu0)
    Nu0 = 1000;
end

if ischar(nToTry)
    nToTry = str2double(nToTry);
end

if ischar(nRequests)
    nRequests = str2double(nRequests);
end

% display(['evaluating: ' num2str(nToTry)]);

X = dlmread(featureFile); %should produce X

if size(X,2)==3 && max(max(X))>1 %sparse format - rejig this to be a matrix

    nDocs = X(1,1);
    nFeat = X(2,1);
    X = X(3:end,:);
    
    % remove noise (arbitrary level)
    threshold = 0;
%     display(['filtering counts below ' num2str(threshold) '; removing ' num2str(sum(X(:,3)<threshold)) ]);
    X( X(:,3)<threshold, 3 ) = 0; 
%     X( X(:,3)>0, 3 ) = 1;

    Xcols{1} = X(X(:,3)>=threshold, 2)+1;
    Xcols{2} = X(X(:,3)>=threshold, 1)+1;
    Xcols{3} = X(X(:,3)>=threshold, 3);
    
    %add one to indexes since matlab indexes from 1
    X = sparse(X(:,1)+1, X(:,2)+1, X(:,3), nDocs,nFeat); 
else
    nFeat = size(X,2);
    nDocs = size(X,1);
    display('VB version needs to convert the matrix - not yet implemented');
end

% T = dlmread(labelFile); %should produce T. For docs with no label, the whole row is zero.
Tcols = dlmread(labelFile);

fid = fopen(topicDocNoPairsFile);
topicDocNoCells = textscan(fid, '%d %s', 'Delimiter', '," ', 'MultipleDelimsAsOne',true);

fid = fopen(fileMapFile);
fileMapCells = textscan(fid, '%d %s %s %s', 'Delimiter', '/_., ');

fileIdxs = fileMapCells{1} + 1; %starts at zero
fileMapNames = fileMapCells{3};

origTsorted = [0; sort(unique(topicDocNoCells{1}))]; 
nClasses = length(origTsorted);

if min(Tcols(:,3))==0 
    %not yet converted Tcols to start from 1
    for d=1:length(Tcols(:,3))   
        Tcols(d,3) = find(origTsorted==Tcols(d,3)); %assuming that the 0 values mean "none of the above"
        if Tcols(d,5) > 0
            Tcols(d,5) = find(origTsorted==Tcols(d,5));
        end
    end    
elseif min(Tcols(:,3))>nClasses
    %not yet converted Tcols to start from 1
    for d=1:length(Tcols(:,3))   
        Tcols(d,3) = find(origTsorted==Tcols(d,3));
        if Tcols(d,5) > 0
            Tcols(d,5) = find(origTsorted==Tcols(d,5));
        end
    end
end

if strcmp(classifierType, 'two-stage')
    confWeight = 1;
elseif nargin < 14    
    %set to one if you don't want to use it
    confWeight = 0.7;
elseif strcmp(confWeight, 'learn')
    
    notConfIdxs = find(Tcols(:,4)==0);    
    confIdxs = find(Tcols(:,4)==1);
    
    nCorrectConf = sparse(1, Tcols(confIdxs,1), double(Tcols(confIdxs,3)==Tcols(confIdxs,5)) );
    nIncConf = sparse(1, Tcols(confIdxs,1), double(Tcols(confIdxs,3)~=Tcols(confIdxs,5)));
    pCorrectConf = (nCorrectConf + 9) ./ (nIncConf + 10 + nCorrectConf);
    
    nCorrectNConf = sparse(1, Tcols(notConfIdxs,1), double(Tcols(notConfIdxs,3)==Tcols(notConfIdxs,5)));
    nIncNConf = sparse(1, Tcols(notConfIdxs,1), double(Tcols(notConfIdxs,3)~=Tcols(notConfIdxs,5)));
    pCorrectNConf = (nCorrectNConf + 5.5) ./ (nIncNConf + 10 + nCorrectNConf );    
    
    confWeight_byWorker = pCorrectNConf ./ pCorrectConf;
    confWeight_byWorker(confWeight_byWorker>1) = 1;
    confWeight = confWeight_byWorker(Tcols(notConfIdxs,1));
end

confWeightedResponses = Tcols(:, 4);
confWeightedResponses(confWeightedResponses==0) = confWeight;

if nargin < 15
    crowdTrust = [20 80];
end

if nargin < 9 || isempty(Alpha0)
    Alpha0 = 1 .* ones(nClasses, 2);
end

labIdxs = Tcols(:,2);
T = zeros(length(labIdxs), nClasses);
T( sub2ind(size(T), (1:length(labIdxs))',Tcols(:,3)) ) = confWeightedResponses;

%old method when we had only one, trusted label per data point -- see below
% unlabIdxs = find(sum(T,2)==0);
% labIdxs = sum(T,2)>0; %indexes of labelled data points

if ~exist('nRequests','var')
    nRequests = 10; %this is currently ignored and fixed anyway
end
if ~exist('nToTry','var')
    nToTry = 50;
end

currentTrain = X(labIdxs,:);
currentLabels = T;
currentTest = X;% previously we trusted the labels and did this: X(unlabIdxs,:);
%we now wish to re-classify all points based on their features -> doesn't
%really use the crowdsourced labels very well -> use VB to fix this

bccSettings = settings.BccSettings();
bccSettings.Alpha = Alpha0;
bccSettings.nu{1} = Nu0;

if strcmp(classifierType, 'VB')
    
    X = X ./ max(max(X,[],1));
    Tvec = zeros(1,nDocs); %ends up just picking one of the classes if there are multiple assignments
    Tvec(Tcols(:,2)) = Tcols(:,3);
    targetWeights = ones(1,nDocs); 
    targetWeights(Tcols(:,2)) = confWeightedResponses;
    
    combiner = combiners.bcc.CIbccVb(bccSettings, nFeat, 1, Tvec, [], nClasses, 2);
    combiner.useLikelihood = true;
    combiner.useNegFeatures = false;
    combiner.targetWeights = targetWeights;
    combiner.setTargets(Tvec);
    [combinedPost Alpha] = combiner.combineDecisions(X');
    P = combinedPost';
    
elseif strcmp(classifierType, 'VB-uncertLabels')
    
    X = X ./ max(max(X,[],1));
    Tvec = sparse(1, nDocs); %have no fully reliable labels except test items
    Tvec(sub2ind(size(Tvec), ones(size(Tcols,1),1), Tcols(:,2))) = Tcols(:,5);
    
    %turn the labels into an additional column for each class. Ones in the 
    %column where the user has selected that class.
    newCol = sparse(Tcols(:,2), Tcols(:,3), confWeightedResponses, size(X,1), nClasses);
    newCol = full(newCol);
    newCol( sum(newCol,2)==0, :) = -1;
    X = [X newCol];
    
    %set an appropriate prior for the crowd's labels
    Alpha0 = repmat(Alpha0, [1, 1, size(X,2)]);
    Alpha0(:,:,nFeat+1:end) = crowdTrust(1);
    for j=1:nClasses
        Alpha0(j, 2, nFeat+j) = crowdTrust(2);
    end
    bccSettings.Alpha = Alpha0;
    combiner = combiners.bcc.CIbccVb(bccSettings, nFeat+nClasses, 1, Tvec, [], nClasses, 2);
    combiner.useLikelihood = true;
    combiner.useNegFeatures = false;    
    [combinedPost Alpha] = combiner.combineDecisions(X');
    P = combinedPost';    
    
elseif strcmp(classifierType, 'VB-workerUncert')
    
    workerMap = unique(Tcols(:,1));
    for w=1:size(Tcols,1)
        Tcols(w,1) = find(Tcols(w,1)==workerMap);
    end
    nWorkers = length(workerMap);
    
    X = X ./ max(max(X,[],1));
    Tvec = sparse(1, nDocs); %have no fully reliable labels except test items
    Tvec(sub2ind(size(Tvec), ones(size(Tcols,1),1), Tcols(:,2))) = Tcols(:,5);
    
    %turn the labels into an additional column for each class for each worker. 
    % Ones in the column where the user has selected that class.
    newCol = sparse(Tcols(:,2), Tcols(:,3)+(Tcols(:,1)-1).*nClasses, confWeightedResponses, size(X,1), nWorkers.*nClasses);
    newCol = full(newCol);
    for w=0:nWorkers-1
        noRespIdxs = sum(newCol(:, (1:nClasses)+(w*nClasses)), 2)==0;
        newCol(noRespIdxs, (w*nClasses)+(1:nClasses)) = -1;
    end
    X = [X newCol];
    
    %set an appropriate prior for the crowd's labels
    Alpha0 = repmat(Alpha0, [1, 1, size(X,2)]);
    Alpha0(:,:,nFeat+1:end) = crowdTrust(1);
    for w=0:nWorkers-1
        for j=1:nClasses
            Alpha0(j, 2, nFeat+j+(w*nClasses)) = crowdTrust(2);
        end    
    end
    bccSettings.Alpha = Alpha0;
    combiner = combiners.bcc.CIbccVb(bccSettings, nFeat+(nClasses*nWorkers), 1, Tvec, [], nClasses, 2);
    combiner.useLikelihood = true;
    combiner.useNegFeatures = false;       
    [combinedPost Alpha] = combiner.combineDecisions(X');
    P = combinedPost';
    
else
    if strcmp(classifierType, 'two-stage-G') || strcmp(classifierType, 'two-stage-G-inclLabels')
        G = 0.9;
        G_workers = trust(labelFile, G);
        G_labels = G_workers(Tcols(:,1));
    else
        G_labels = ones(size(Tcols(:,1)));
    end
    if nargin < 9
        [P_class,P_feat_class,Alpha] = multi_classifier(G_labels, currentTrain,currentLabels, Alpha0);
    else
        [P_class,P_feat_class,Alpha] = multi_classifier(G_labels, currentTrain,currentLabels, Alpha0, Nu0);
    end
    P = multi_classify(P_class,P_feat_class,currentTest,true,false);
    
    if strcmp(classifierType,  'two-stage-inclLabels') || strcmp(classifierType, 'two-stage-G-inclLabels')
        pLabCorrect = 0.7; 
        pIncorrect = 1-pLabCorrect;
        
        T = sparse( Tcols(:,2), Tcols(:,3), confWeightedResponses, size(P,1), size(P,2) );
        labelCounts = repmat(sum(T,2), 1, size(T,2));
        T(T==0) = labelCounts(T==0) .* log(pIncorrect);
        T(T>0) = labelCounts(T>0) .* log(pLabCorrect);
        
        P = exp(log(P) + T); 
        P = P ./ repmat( sum(P,2), 1, size(P,2) ); %incorporate labels
    end
end

testResults = cell(size(P,2),1);
for j=1:nClasses
    testResults{j} = P(:, j)';
end
if strcmp(selectMethod, 'uncertGreedy')
    samplesToRequest = labelSelectUncertGreedy(nRequests, X, P, Alpha, Tcols, P_feat_class, P_class);    
elseif strcmp(selectMethod, 'uncertClust')
    samplesToRequest = labelSelectUncertDiff(nToTry, nRequests, X, P, Alpha, Tcols);
elseif strcmp(selectMethod, 'random')
    samplesToRequest = randi(size(X,1), nRequests, 1);
elseif strcmp(selectMethod, 'uncert')
    samplesToRequest = labelSelectUncert(nRequests, X, P, Alpha, Tcols);
end
if exist('samplesToRequest', 'var')
    dlmwrite(outputRequestFile, samplesToRequest);
else
    samplesToRequest = [];
end

%OUTPUTTING IN THE TREC SUBMISSION FORMAT

%choose thresholds from ROC
binaryRes = zeros(size(P));

threshold = zeros(1,nClasses);
labIdxs = unique(labIdxs);
T = sparse(Tcols(:,2), Tcols(:,3), 1, size(P,1), nClasses);
for j=2:size(P,2)     
    [tpr fpr theta] = roc(T(labIdxs,j)', P(labIdxs,j)');
    rating = (tpr + (1-fpr)) .* (tpr>fpr);
    [sortRatings sortIdx] = sort(rating, 'ascend');
    maxRatIdx = sortIdx(end);
%     [maxRating maxRatIdx] = max(rating);
    if theta(maxRatIdx) < 1/nClasses
        threshold(j) = (theta(maxRatIdx) + 1./nClasses) / 2;
    else
        threshold(j) = theta(maxRatIdx);
    end
    if threshold(j)==1
        threshold(j) = 1./nClasses;
    end
end    
threshMat = repmat(threshold, size(P,1), 1);
maxPs = repmat( max(P,[],2), 1, size(P,2));
binaryRes(P>threshMat) = 1;% &  P>=maxPs) = 1;
    
for j=1:nClasses
    if sum(binaryRes(:,j),1) > nDocs * 0.1 || sum(binaryRes(:,j),1) < nDocs * 0.01
        [vals idxs] = sort(P(:,j), 'ascend');
        threshold(j) = vals(round(nDocs*0.9));

        threshMat = repmat(threshold, size(P,1), 1);
        maxPs = repmat( max(P,[],2), 1, size(P,2));
        binaryRes(:,j) = 0;
        binaryRes(P>=maxPs) = 1;%P(:,j)>threshold(j), j) = 1;% &  
    end
end
display(['Binary decision threshold for each class: ' num2str(threshold)]);

pairsRes = zeros(size(topicDocNoCells{1}, 1),4);

if ~writeOutput
    return
end

fid = fopen(sprintf(outputResultFile,runName), 'w');
for line=1:size(topicDocNoCells{1})
    
    topic = topicDocNoCells{1}(line);
    topicIdx = find(origTsorted==topic);
    
    docNo = topicDocNoCells{2}(line); docNo = docNo{1};
    docIdx = fileIdxs(( strcmp(docNo, fileMapNames) ));
        
    if isempty(docIdx) 
        pVal = P_class(topicIdx);
        binVal = 0;
        display(['Could not find the document: ' docNo]);
    elseif topicIdx>size(P,2)
        display(['TopicIdx invalid: ' num2str(topicIdx) ', max topicIdx is ' num2str(size(P,2)) ]);
        binVal = 0;
        pVal = 0.2;
    else
        binVal = binaryRes(docIdx, topicIdx);  
        if sum(binVal)>=1
            binVal = 1;
        else
            binVal = 0;
        end
                
        pVal = P(docIdx(1), topicIdx);     
%         [maxP, maxIdx] = max(P(docIdx,:), [], 2);
%         if maxIdx==topicIdx
%             binVal = 1;
%         else
%             binVal = 0;
%         end
    end
    outputString = sprintf('%d %s %d %4f %s\n', topic, docNo, binVal, pVal, runName);
    pairsRes(line,3) = binVal;
    pairsRes(line,4) = pVal;
    pairsRes(line,1) = topic;
    pairsRes(line,2) = docIdx;
    fprintf(fid, outputString);
end
fclose(fid);


