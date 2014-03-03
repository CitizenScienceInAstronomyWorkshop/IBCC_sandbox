%Select labels using the TrecLabelRequesting method (rough approximation to
%selecting combinations by IG, uses Pseudo count method to evaluate) 

%using pseduo-count update method with a separate training phase instead of
%VB. VB could be used later once labels have been selected.

% featureFile = '/homes/49/edwin/data/trec/antonio/matrices/matrices_DFR.txt';
% featureFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-2000-Matrix/matrices_2000_topic.txt';
% fileMapFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-2000-Matrix/fileMap.txt';
% featureFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-1500-Matrix/matrices_1500_topic.txt';
% fileMapFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-1500-Matrix/fileMap.txt';

X = dlmread(featureFile, ' ');
nDocs = X(1,1);
nFeatures = X(2,1);

% fid = fopen('/homes/49/edwin/data/trec/antonio/matrices/fileMap.txt');


fid = fopen(fileMapFile);
% fileMapCells = textscan(fid, '%d%s%s%s%s%s%s%s%s%s%s %s %d %d %s', 'Delimiter', '/_. ');
fileMapCells = textscan(fid, '%d%s %s %d %d %s', 'Delimiter', '/_. ');

fileIdxs = fileMapCells{1} + 1; %starts at zero
fileMapTopics = fileMapCells{4};%fileMapCells{13};
fileMapRel = fileMapCells{5};%fileMapCells{14};

% fileMapTopics = fileMapCells{13};
% fileMapRel = fileMapCells{14};

Tcols = zeros(0,3);

nQrels = length(fileMapRel);
for d=1:nQrels
    if fileMapRel(d) == 1
        Tcols = [Tcols; [0 fileIdxs(d) fileMapTopics(d)]];
    end
end

nClasses = length(unique(Tcols(:,3)))+1;
Alpha0 = 1 .* ones(nClasses, 2);
Nu0 = 1;

Tsorted = sort(unique(Tcols(:,3)));

for d=1:length(Tcols(:,3))   
    Tcols(d,3) = find(Tsorted==Tcols(d,3))+1;
end

T = sparse(double(Tcols(:,2)), double(Tcols(:,3)), 1, nDocs, length(Tsorted)+1);

knownIdxs = find(sum(T,2)>0);

%treat unknown ones as if they are definitely none of the above
unknownIdxs = find(sum(T,2)==0);
T(unknownIdxs,1) = 1;
Tcols = [Tcols; [zeros(length(unknownIdxs),1) unknownIdxs ones(length(unknownIdxs),1)]];

nRounds = 10;
aucs = zeros(nClasses, nRounds);
lam_uncert = zeros(1,nRounds);

nRequests = 50;
nToTry = 250;

labelFile = '/homes/49/edwin/data/trec/trec7Matrices/labelFile';
outputRequestFile = '/homes/49/edwin/data/trec/trec7Matrices/outputRequestFile';
outputResultFile = '/homes/49/edwin/data/trec/trec7Matrices/outputResultFile.csv';
topicDocNoPairsFile = '/homes/49/edwin/data/trec/qrels/topic_docno_pairs_trec7_relevantOnly.csv';

%write currentLabels to labelFile

if ~exist('initLabelIds', 'var')
    initLabelIds = randperm(size(Tcols,1)); 
    initLabelIds = initLabelIds(1:500);

    currentLabels = Tcols(initLabelIds,:);    
    
    %corrupt labels to simulate uncertain decision makers
    pIncorrect = 0.2;
    corruptLabels = find(rand(1,length(initLabelIds)) < pIncorrect);
    changeToLabels = randi(nClasses-1, length(corruptLabels), 1);
    currentLabels(corruptLabels,3) = double(currentLabels(corruptLabels,3)) + changeToLabels;
    currentLabels(currentLabels(:,3)>nClasses,3) = currentLabels(currentLabels(:,3)>nClasses,3) - nClasses;
end

for f=1:nRounds   
    dlmwrite(labelFile, currentLabels);
    
    display(num2str(f));
  
    
    [samplesToRequest, P, binaryRes] = classifyAndSelectDoc_greedy( nToTry, nRequests, ...
        featureFile, labelFile, outputRequestFile, outputResultFile, ...
        topicDocNoPairsFile, fileMapFile, [], [], false);

    pTrue = cell(1, nClasses);
    binLabs = cell(1, nClasses);
        
    [maxCert maxClass] = max(P,[],2);
%     binaryRes = zeros(size(P));
      
%     threshold = zeros(1,nClasses);
    
%     for j=2:size(P,2)     
%         [tpr fpr theta] = roc(T(:,j)', P(:,j)');
%         rating = tpr + (1-fpr) .* (tpr>fpr);
%         [maxRating maxRatIdx] = max(rating);
%         threshold(j) = theta(maxRatIdx);
%     end    
%     threshMat = repmat(threshold, size(P,1), 1);
%     binaryRes( P>threshMat ) = 1;
%     binaryRes(sub2ind(size(P),1:size(P,1),maxClass')) = 1;
%     binaryRes(repmat(maxCert, 1, nClasses)-P<0.07 & P>(1/nClasses)) = 1;
    
    nFp = 0.5; nPos = 1;
    nFn = 0.5; nNeg = 1;
    
    for d=1:nQrels
        
        topic = find(Tsorted==fileMapTopics(d))+1;
        doc = fileIdxs(d);
        
        pTrue{topic} = [pTrue{topic} P(doc,topic)];
        binLabs{topic} = [binLabs{topic} fileMapRel(d)];        
        
        if fileMapRel(d)==1
                       
            nPos = nPos + 1;
            
%             if P(doc, topic) < 0.2
            if binaryRes(doc,topic) ~= 1
                nFn = nFn + 1;
            end
        else
            nNeg = nNeg + 1;
            
%             if P(doc, topic) >= 0.2
            if binaryRes(doc,topic) == 1
                nFp = nFp + 1;
            end
        end
    end
    
    testResults = cell(size(P,2),1);
    for j=2:size(P,2)        
        auc = graphs.ClassifierPerformanceGraph.drawRoc(pTrue{j}, ...
            binLabs{j}, {'pseudo ibcc'}, false, true);
        aucs(j,f) = auc;
    end
        
    fpr = nFp ./ nNeg;
    fnr = nFn ./ nPos;
    
    lam_uncert_f = (log(fpr./(1-fpr)) + log(fnr./(1-fnr)))./2;
    lam_uncert_f = exp(lam_uncert_f) ./ (1+exp(lam_uncert_f));
    lam_uncert(f) = lam_uncert_f;

    if f<nRounds
        tmp = Tcols(ismember(Tcols(:,2),samplesToRequest), :);
        currentLabels = [currentLabels; tmp];
    end
end

lam_uncert(f)
aucs(:,f)'