function [ aucs ] = evaluateResults( P, binaryRes, Tsorted, fileMapFile, topicDocNoPairsFile, qRels, noGraphs )
%EVALUATERESULTS Summary of this function goes here
%   Detailed explanation goes here

    nClasses = size(P,2);

    %For the test pairs from each topic, probability of "true" match
    %according to our algorithm
    pTrue = cell(1, nClasses);
    %Ground truth binary labels for the same pairs
    binLabs = cell(1, nClasses);
            
    nFp = 0.5; nPos = 1;
    nFn = 0.5; nNeg = 1;
        
    fid = fopen(fileMapFile);
    fileMapCells = textscan(fid, '%d %s %s %s', 'Delimiter', '/_. ');
    fileIdxs = fileMapCells{1} + 1; %starts at zero
    filenames = fileMapCells{3};

    fid = fopen(topicDocNoPairsFile);
    topicDocNoCells = textscan(fid, '%d %s', 'Delimiter', '," ', 'MultipleDelimsAsOne',true);
    pairsTopics = topicDocNoCells{1};
    pairsFiles = topicDocNoCells{2};
    pairsFileIdxs = zeros(size(pairsFiles)); 
    nPairs = length(pairsTopics);
    
    %Produce arrays of our results and ground truth results -------
    for d=1:nPairs
        pairsToFileMapIdx = strcmp(pairsFiles(d),filenames);
        if sum(pairsToFileMapIdx)==0
            continue;
        end
        pairsFileIdxs(d) = fileIdxs(pairsToFileMapIdx);

        topic = find(Tsorted==pairsTopics(d));
        doc = pairsFileIdxs(d);
        
        pTrue{topic} = [pTrue{topic} P(doc,topic)];
        binLabs{topic} = [binLabs{topic} qRels(doc,topic)];        
%         
%         if qRels(doc,topic)==1         
%             nPos = nPos + 1;
%             if binaryRes(doc,topic) ~= 1
%                 nFn = nFn + 1;
%             end
%         else
%             nNeg = nNeg + 1;
%             if binaryRes(doc,topic) == 1
%                 nFp = nFp + 1;
%             end
%         end
    end
    
    %AUC curves ---------------------------------------------------
    aucs = zeros(size(P,2)-1,1);
    for j=2:size(P,2)        
        auc = graphs.ClassifierPerformanceGraph.drawRoc(pTrue{j}, ...
            binLabs{j}, {'topic'}, false, noGraphs, false);
        aucs(j-1) = auc;
    end
        
    %LAM ---------------------------------------------------------
%     fpr = nFp ./ nNeg;
%     fnr = nFn ./ nPos;
%     
%     lam = (log(fpr./(1-fpr)) + log(fnr./(1-fnr)))./2;
%     lam = exp(lam) ./ (1+exp(lam));
end

