function [qrels aucs lam] = evaluateTrecResults(pairsRes, fileMapFile, qrels, combineClasses)

%The option combineClasses is a flag to indicate whether the ROC should be
%calculated for all classes combined (as in the challenge itself) or for
%each separately.

    if nargin < 4
        combineClasses = false;
    end

    if nargin < 3 || isempty(qrels)
        qrelFile = '/homes/49/edwin/data/trec/qrels/qrels.trec8.adhoc.parts1-5';
        qrelCols = textscan(fopen(qrelFile), '%d%d%s%d', 'Delimiter', ' ');
        definiteNegs = qrelCols{4}==0;
        qrelCols{4}(definiteNegs) = -1;
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
        qrelTopics = qrelCols{1}(validIdxs);
        qrelDocs = qrelDocIdxs(validIdxs);
        qrelRel = qrelCols{4}(validIdxs);

        qrels = sparse(double(qrelTopics), qrelDocs, double(qrelRel)); %topics, documents -> binary 1 or -1; 0s means no entry
    end
    topicMap = unique(pairsRes(:,1));
    nClasses = length(topicMap);
    
    nPairs = size(pairsRes,1);
    
    pTrue = cell(1, nClasses);
    binLabs = cell(1, nClasses);
            
    nFp = 0.5; nPos = 1;
    nFn = 0.5; nNeg = 1;
   
    nSkipped = 0;
    
    for d=1:nPairs%size(binaryRes,1)
%         for t=1:(size(binaryRes,2)-1)        
            realTopicId = pairsRes(d,1);
            topicIdx = find( realTopicId==topicMap );
            doc = pairsRes(d,2);
%             topic = find(t+1==topicMap);

            if qrels(realTopicId,doc)==0
                display('no true value for pair');
                nSkipped = nSkipped + 1;
                continue
            elseif qrels(realTopicId,doc)==-1
                trueLab = 0;
            else
                trueLab = 1;
            end
            if combineClasses
                pTrue{1} = [pTrue{1} pairsRes(d,4)];        
                binLabs{1} = [binLabs{1} trueLab];
            else
                pTrue{topicIdx} = [pTrue{topicIdx} pairsRes(d,4)];        
                binLabs{topicIdx} = [binLabs{topicIdx} trueLab];
            end

            if qrels(realTopicId, doc)==1         
                nPos = nPos + 1;
                if pairsRes(d,3) ~= 1
                    nFn = nFn + 1;
                end
            else
                nNeg = nNeg + 1;
                if pairsRes(d,3) == 1
                    nFp = nFp + 1;
                end
            end
%         end
    end
    
    display(['number of pairs with no true label: ' num2str(nSkipped)]);
    
    aucs = zeros(1, nClasses);    
    figure;
    
    if combineClasses
        nClasses = 1;
    end
    
    for j=1:nClasses    
        auc = graphs.ClassifierPerformanceGraph.drawRoc(pTrue{j}, ...
            binLabs{j}, {'ibcc'}, false, false, false);
        aucs(j) = auc;
    end
        
    fpr = nFp ./ nNeg;
    fnr = nFn ./ nPos;
    
    lam = (log(fpr./(1-fpr)) + log(fnr./(1-fnr)))./2;
    lam = exp(lam) ./ (1+exp(lam));
    
    title('TREC 8 Results');
