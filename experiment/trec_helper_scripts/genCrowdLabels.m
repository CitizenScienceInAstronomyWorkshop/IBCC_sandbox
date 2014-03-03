function [crowdLabels, currentLabels, qRels] = genCrowdLabels( fileMapFile, ...
    topicDocNoPairsFile, qRels, nInitialLabels, nWorkers, nClasses )
%GENCROWDLABELS Generate labels for the crowd by corrupting qrels or making
%the qrels up if they are missing.
    
    %ORACLE LABELS FOR ALL KNOWN PAIRS -----------------------------------------
    %For anything we have ground truth for, we can ask for a corrupted version
    %of that.
    %Load File map
    fid = fopen(fileMapFile);
    fileMapCells = textscan(fid, '%d %s %s %s', 'Delimiter', '/_. ');
    fileIdxs = fileMapCells{1} + 1; %starts at zero

    fid = fopen(topicDocNoPairsFile);
    topicDocNoCells = textscan(fid, '%d %s', 'Delimiter', '," ', 'MultipleDelimsAsOne',true);
    pairsTopics = topicDocNoCells{1};

    if ~exist('qRels', 'var')
        nTopics = length(unique(pairsTopics)) + 1; % +1 for the none-of-the-above category
        nDocs = length(unique(fileIdxs));
        
        qRels = randi(2, nDocs, nTopics) - 1; %Got to get these from the file if we want to test
    else
        nDocs = size(qRels, 1);
    end    
    
    Tcols = zeros(nDocs,5);
    Tcols(:,2) = (1:nDocs)';
    %treat unknown ones as if they are definitely none of the above
    [~, trueLabels] = max(qRels, [], 2);
    Tcols(:,3) = trueLabels;
    % Tcols(:,5) = trueLabels;

    %CROWD LABELS FROM ORACLE LABELS ------------------------------------------
    %write currentLabels to labelFile
    nKnownPairs = size(Tcols,1);

    %Set of label IDs we know before requesting any new ones
    initLabelIds = randperm(nKnownPairs); 
    initLabelIds = initLabelIds(1:nInitialLabels);
    %initLabelIds = []; % known nothing from crowd in first round - set inf. priors

    crowdLabels = Tcols;

    %1: Set Worker IDs
    crowdLabels(:,1) = randi(nWorkers, nKnownPairs, 1);

    %2: Doc IDs - don't change

    %3: Set workers' responses
    %corrupt labels to simulate uncertain decision makers
    pIncorrect = rand(nWorkers, 1) .* 0.3;
    pIncorrect = pIncorrect(crowdLabels(:,1));
    corruptLabels = find(rand(length(crowdLabels),1)-pIncorrect < 0);
    changeToLabels = randi(nClasses-1, length(corruptLabels), 1);
    crowdLabels(corruptLabels,3) = double(crowdLabels(corruptLabels,3)) + changeToLabels;
    crowdLabels(crowdLabels(:,3)>nClasses,3) = crowdLabels(crowdLabels(:,3)>nClasses,3) - nClasses;

    %4: set confidence labels randomly
    crowdLabels(:,4) = randi(2, nKnownPairs, 1) - 1;

    %5: true labels - don't change

    currentLabels = crowdLabels(initLabelIds,:);
end

