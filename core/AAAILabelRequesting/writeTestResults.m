function [binaryRes pairsRes] = writeTestResults( P, combiner, fileMapFile, outputResultFile, ...
    topicDocNoPairsFile, writeOutput, runName, binaryMethod, Tcols, origTsorted )
%WRITETESTRESULTS Write out the results for the TREC test pairs
%(documents/topic pairs) and the binary results.
%   Binary and pairs results - write to file in TREC submission format.

if ~exist('runName', 'var')
    runName = 'Orc2013';
end

if ~exist('writeOutputs', 'var')
    writeOutput = true;
end

%binary method can be set to choose the binary labels from a set of
%probabilities in different ways. 
% 'posterior' - Choose the most probable class for a document
% 'likelihood' - Choose the most likely (ignoring class priors) for a doc
% 'proportions' - Choose the overall number of positives for each class to match the class priors
% 'roc' - Used in the TREC submission to try to find the correct threshold
% for each class. 
if ~exist('binaryMethod', 'var')
    binaryMethod = 'proportions';
end

fid = fopen(topicDocNoPairsFile);
topicDocNoCells = textscan(fid, '%d %s', 'Delimiter', '," ', 'MultipleDelimsAsOne',true);

%%%%%%%%%%%% Choosing Binary Labels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nDocs = size(P,1);

%If this is time-consuming, only needs to be run at the end before
%submission
%Can now be replaced/augmented with some evaluation metrics using final
%released judgements

binaryResSet = {};
binaryRes = zeros(size(P));

if strfind(binaryMethod, 'proportions')
    Kappa = combiner.Nu ./ sum(combiner.Nu);
    nPerClass = round(Kappa .* nDocs);
    for j=2:size(P,2)     
        [~, sortedDocIdxs] = sort(P(:,j), 'descend');
        chosen = sortedDocIdxs(1:nPerClass(j));
        binaryRes(chosen, j) = 1;
    end
    binaryResSet = [binaryResSet binaryRes];
end
if strfind(binaryMethod, 'posterior')
    [~, chosenClasses] = max(P, [], 2);
    binaryRes( sub2ind(size(binaryRes), 1:nDocs, chosenClasses') );
    binaryResSet = [binaryResSet binaryRes];    
end
if strfind(binaryMethod, 'likelihood')
    Kappa = repmat(combiner.Nu ./ sum(combiner.Nu), nDocs, 1);
    likelihood = P ./ Kappa;
    [~, chosenClasses] = max(likelihood, [], 2);
    binaryRes( sub2ind(size(binaryRes), 1:nDocs, chosenClasses') );    
    binaryResSet = [binaryResSet binaryRes];
end
if strfind(binaryMethod, 'roc')
    nClasses = size(P,2);
    threshold = zeros(1,nClasses);
    labIdxs = unique(Tcols(:,2));    
    T = sparse(Tcols(:,2), Tcols(:,3), 1, nDocs, nClasses);
    
    for j=2:size(P,2)     
        [tpr fpr theta] = roc(T(labIdxs,j)', P(labIdxs,j)');
        rating = (tpr + (1-fpr)) .* (tpr>fpr);
        [~, sortIdx] = sort(rating, 'ascend');
        maxRatIdx = sortIdx(end);
        if theta(maxRatIdx) < 1/nClasses
            threshold(j) = (theta(maxRatIdx) + 1./nClasses) / 2;
        else
            threshold(j) = theta(maxRatIdx);
        end
        if threshold(j)==1
            threshold(j) = 1./nClasses;
        end
    end    
    threshMat = repmat(threshold, nDocs, 1);
    binaryRes(P>threshMat) = 1;

    for j=1:nClasses
        if sum(binaryRes(:,j),1) > nDocs * 0.1 || sum(binaryRes(:,j),1) < nDocs * 0.01
            vals = sort(P(:,j), 'ascend');
            threshold(j) = vals(round(nDocs*0.9));

            maxPs = repmat( max(P,[],2), 1, size(P,2));
            binaryRes(:,j) = 0;
            binaryRes(P>=maxPs) = 1;
        end
    end
%     display(['Binary decision threshold for each class: ' num2str(threshold)]);    
    binaryResSet = [binaryResSet binaryRes];
end

if ~writeOutput
    return
end

%%%%%%%%% Writing the results for the test pairs %%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(fileMapFile);
fileMapCells = textscan(fid, '%d %s %s %s', 'Delimiter', '/_., ');

fileIdxs = fileMapCells{1} + 1; %starts at zero
fileMapNames = fileMapCells{3};

pairsRes = zeros(size(topicDocNoCells{1}, 1),4);

fid = fopen(sprintf(outputResultFile,runName), 'w');
for line=1:size(topicDocNoCells{1})
    
    topic = topicDocNoCells{1}(line);
    topicIdx = find(origTsorted==topic);
    
    docNo = topicDocNoCells{2}(line); docNo = docNo{1};
    docIdx = fileIdxs(( strcmp(docNo, fileMapNames) ));
        
    if isempty(docIdx) 
%         pVal = Kappa(topicIdx);
%         binVal = 0;
%         display(['Could not find the document: ' docNo]);
        continue;
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
    end
    outputString = sprintf('%d %s %d %4f %s\n', topic, docNo, binVal, pVal, runName);
    fprintf(fid, outputString);

    pairsRes(line,1) = topic;
    pairsRes(line,2) = docIdx(1);
    pairsRes(line,3) = binVal;
    pairsRes(line,4) = pVal;
end
fclose(fid);

