%Select labels using the TrecLabelRequesting method (rough approximation to
%selecting combinations by IG, uses Pseudo count method to evaluate) 

%using pseduo-count update method with a separate training phase instead of
%VB. VB could be used later once labels have been selected.

% featureFile = '/homes/49/edwin/data/trec/antonio/matrices/matrices_DFR.txt';
% featureFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-2000-Matrix/matrices_2000_topic.txt';
% fileMapFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-2000-Matrix/fileMap.txt';
% featureFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-1500-Matrix/matrices_1500_topic.txt';
% fileMapFile = '/homes/49/edwin/data/trec/antonio/TopicCorpus-1500-Matrix/fileMap.txt';

if ~exist('selectMethod', 'var')
    selectMethod = 'uncertClust';
end

if ~exist('classifierMethod', 'var')
    classifierMethod = 'two-stage';
end

if ~exist('crowdTrust', 'var')
    crowdTrust = [20 80];
end

labelFile = '/homes/49/edwin/data/trec/testOutput/labelFile_selmethods';
outputRequestFile = '/homes/49/edwin/data/trec/testOutput/outputRequestFile';
outputResultFile = '/homes/49/edwin/data/trec/testOutput/outputResultFile.csv';

topicDocNoPairsFile = '/homes/49/edwin/data/trec/qrels/topic_docno_pairs_trec7_relevantOnly.csv';
% topicDocNoPairsFile =
% '/homes/49/edwin/data/trec/trec8_sample/trec2012-crowdsourcing-text-task.topic-docnos.txt';
% %use this file to test whether we produce real output. The topics are
% different though so will all have 0s against them.

X = dlmread(featureFile, ' ');
nDocs = X(1,1);
nFeatures = X(2,1);

fid = fopen(fileMapFile);
% fileMapCells = textscan(fid, '%d%s%s%s%s%s%s%s%s%s%s %s %d %d %s', 'Delimiter', '/_. ');
fileMapCells = textscan(fid, '%d%s %s %d %d %s', 'Delimiter', '/_. ');

fileIdxs = fileMapCells{1} + 1; %starts at zero
fileMapTopics = fileMapCells{4};%fileMapCells{13};
fileMapRel = fileMapCells{5};%fileMapCells{14};


Tcols = zeros(0,5);

nQrels = length(fileMapRel);
for d=1:nQrels
    if fileMapRel(d) == 1
        Tcols = [Tcols; [0 fileIdxs(d) fileMapTopics(d) randi(2)-1 fileMapTopics(d)]]; %generate a random confidence label just to check the code works
    end
end

nClasses = length(unique(Tcols(:,3)))+1;
Alpha0 = 1 .* ones(nClasses, 2);
Nu0 = 1;

Tsorted = sort(unique(Tcols(:,3)));

for d=1:length(Tcols(:,3))   
    Tcols(d,3) = find(Tsorted==Tcols(d,3))+1;
    Tcols(d,5) = find(Tsorted==Tcols(d,5))+1;
end

T = sparse(double(Tcols(:,2)), double(Tcols(:,3)), 1, nDocs, length(Tsorted)+1);

knownIdxs = find(sum(T,2)>0);

%treat unknown ones as if they are definitely none of the above
unknownIdxs = find(sum(T,2)==0);
T(unknownIdxs,1) = 1;
Tcols = [Tcols; [zeros(length(unknownIdxs),1) unknownIdxs ...
    ones(length(unknownIdxs),1) randi(2,length(unknownIdxs),1)-1 ones(length(unknownIdxs),1)]];

nRounds = 500;
aucs = zeros(nClasses, nRounds);
lam_uncert = zeros(1,nRounds);

nRequests = 10;
nToTry = 200;

%write currentLabels to labelFile

if ~exist('initLabelIds', 'var')
    initLabelIds = randperm(size(Tcols,1)); 
    initLabelIds = initLabelIds(1:1);
   
    Tcols(:,1) = randi(30, size(Tcols,1), 1);

    %corrupt labels to simulate uncertain decision makers
    pIncorrect = rand(30, 1) .* 0.3;
    pIncorrect = pIncorrect(Tcols(:,1));
    
    crowdLabels = Tcols;
    corruptLabels = find(rand(length(crowdLabels),1)-pIncorrect < 0);
    changeToLabels = randi(nClasses-1, length(corruptLabels), 1);
    crowdLabels(corruptLabels,3) = double(crowdLabels(corruptLabels,3)) + changeToLabels;
    crowdLabels(crowdLabels(:,3)>nClasses,3) = crowdLabels(crowdLabels(:,3)>nClasses,3) - nClasses;
    
    currentLabels = crowdLabels(initLabelIds,:);
    initialLabels = currentLabels;
%     currentLabels(1:5000,5) = 0;
else
    currentLabels = initialLabels;
end

% crowdLabels = '/homes/49/edwin/data/trec/trec8_sample/labels.csv';

for f=1:nRounds   
    dlmwrite(labelFile, currentLabels);
    
    display(num2str(f));
  
    [samplesToRequest, P, binaryRes] = classifyAndSelectDoc( nToTry, nRequests, ...
        featureFile, labelFile, outputRequestFile, outputResultFile, ...
        topicDocNoPairsFile, fileMapFile, [], [], false, selectMethod, classifierMethod, 1, crowdTrust);%, Alpha0, Nu0 );

    pTrue = cell(1, nClasses);
    binLabs = cell(1, nClasses);
        
    [maxCert maxClass] = max(P,[],2);
    
    nFp = 0.5; nPos = 1;
    nFn = 0.5; nNeg = 1;
    
    for d=1:nQrels
        
        topic = find(Tsorted==fileMapTopics(d))+1;
        doc = fileIdxs(d);
        
        pTrue{topic} = [pTrue{topic} P(doc,topic)];
        binLabs{topic} = [binLabs{topic} fileMapRel(d)];        
        
        if fileMapRel(d)==1         
            nPos = nPos + 1;
            if binaryRes(doc,topic) ~= 1
                nFn = nFn + 1;
            end
        else
            nNeg = nNeg + 1;
            if binaryRes(doc,topic) == 1
                nFp = nFp + 1;
            end
        end
    end
    
    testResults = cell(size(P,2),1);
    if f==nRounds
        figure;
    end
    for j=2:size(P,2)        
        if f==nRounds
            auc = graphs.ClassifierPerformanceGraph.drawRoc(pTrue{j}, ...
                binLabs{j}, {'topic'}, false, false, false);            
        else
            auc = graphs.ClassifierPerformanceGraph.drawRoc(pTrue{j}, ...
                binLabs{j}, {'topic'}, false, true);
        end
        aucs(j,f) = auc;
    end
        
    fpr = nFp ./ nNeg;
    fnr = nFn ./ nPos;
    
    lam_uncert_f = (log(fpr./(1-fpr)) + log(fnr./(1-fnr)))./2;
    lam_uncert_f = exp(lam_uncert_f) ./ (1+exp(lam_uncert_f));
    lam_uncert(f) = lam_uncert_f;

    if f<nRounds
        tmp = crowdLabels(ismember(Tcols(:,2),samplesToRequest), :);
        currentLabels = [currentLabels; tmp];
    end
end

lam_uncert(f)
aucs(:,f)'