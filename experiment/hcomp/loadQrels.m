function [qrels, qrelsNR] = loadQrels(fileMapFile, qrelFile, topicDocNoPairsFile)

%The option combineClasses is a flag to indicate whether the ROC should be
%calculated for all classes combined (as in the challenge itself) or for
%each separately.

fid = fopen(topicDocNoPairsFile);
topicDocNoCells = textscan(fid, '%d %s', 'Delimiter', '," ', 'MultipleDelimsAsOne',true);
origTsorted = [0; sort(unique(topicDocNoCells{1}))]; 

%columns: 1= topic, 2 = ?, 3=filename, 4=true/false
qrelCols = textscan(fopen(qrelFile), '%d%d%s%d', 'Delimiter', ' ');
% definiteNegs = qrelCols{4}==0;
% qrelCols{4}(definiteNegs) = -1;
qrelDocIdxs = zeros(length(qrelCols{4}),1);    

fid = fopen(fileMapFile);
fileMapCells = textscan(fid, '%d %s %s %s', 'Delimiter', '/_., ');

fileIdxs = fileMapCells{1} + 1; %starts at zero
fileMapNames = fileMapCells{3};

for q=1:length(qrelCols{3})
    docNo = qrelCols{3}(q);
    docIdx = fileIdxs(( strcmp(docNo, fileMapNames) ));
    if isempty(docIdx)
        continue
    end
    qrelDocIdxs(q) = double(docIdx(1));
end

validIdxs = qrelDocIdxs~=0;
qrelDocs = qrelDocIdxs(validIdxs);
qrelTopics = qrelCols{1}(validIdxs);
for t=1:length(qrelTopics)
    topic = find(origTsorted==qrelTopics(t));
    if ~isempty(topic)
        qrelTopics(t) = topic;
    end
end
qrelRel = qrelCols{4}(validIdxs);
qrelRelNR = qrelRel;
qrelRelNR(qrelRel==0) = -1; %confirmed not relevant

%0 corresponds to none-of-the-above category, so this is assumed if there is no entry
qrels = sparse(qrelDocs, double(qrelTopics), double(qrelRel)); 
%contains -1 for the confirmed not relevant topics
qrelsNR = sparse(qrelDocs, double(qrelTopics), double(qrelRelNR)); 
