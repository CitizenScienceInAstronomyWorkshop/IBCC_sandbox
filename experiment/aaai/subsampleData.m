if ~exist('chosenIdx','var')
    chosenIdx = [2 5 6]';%[3 7 9]; %In order of most examples: 8, 10, 4, 9, 3, 5, 6, 2, 7, 11
else
    display(['Using topic indexes: ' num2str(chosenIdx')]);
end

if ~exist('batchSize','var')
    batchSize = 500;
end

for repeat=1:nRepeats
    
    display(['Generating document set for test no. ' num2str(repeat)]);
    
    %File listing the subsample of documents we are actually testing on
    selectedDocsFile = [runRootDir '/selectedDocs_' num2str(repeat) '.mat'];
    
    bootstrapFile = [runRootDir '/bootstrap_CrowdLabels' num2str(repeat) '.mat'];
    
    if exist(selectedDocsFile,'file') && exist(bootstrapFile,'file')
        continue;
    end

    %select topic 4, which has 75 positive examples. Or topic 3 with 45 examples.
    %Increase dataset size? If it doesn't take much longer + intelligent tasking still works on small no. clusters

    %Ideally we do each topic at a time in separate runs
    %Otherwise we're doing multi-class but only testing pairs of
    %"none-of-the-above" versus "confirmed topic".
    %However, this multi-class hack actually makes sense for the TREC search
    %application. We can also show improvements without comparing all possible
    %pairs. This would miss any confusion between topics; in search
    %application it's not necessary to remove this confusion as documents can
    %be in multiple classes.

    nClasses = length(chosenIdx)+1;

    posEgs = [];
    for t=1:length(chosenIdx)
        notChosenIdx = chosenIdx(chosenIdx~=chosenIdx(t));
        newIdxs = find( qRels(:,chosenIdx(t)) & qRels_old(:,chosenIdx(t)) & ...
            sum(qRels(:,notChosenIdx),2)<1 & ... %not really confirmed as negative but don't have enough data for that
            sum(qRels_old(:,notChosenIdx),2)<1 );

        newIdxs(ismember(newIdxs,posEgs)) = [];
        posEgs = [posEgs; newIdxs];
    end

    %select a further 9025 random documents; problem - these are not confirmed
    %as irrelevant. Check with QRELS from TREC8.
    negEgs = find( sum(qRels(:,chosenIdx),2)<1 & ...
        sum(qRels_old(:,chosenIdx), 2)<1 );

    negEgs( ismember(negEgs, posEgs) ) = [];

    posBatchSize = length(posEgs);
    negBatchSize = batchSize - posBatchSize;
    startIdx = (repeat-1) * negBatchSize + 1;

    nMissing = negBatchSize;
    negEgsSub = [];

    while nMissing > 0
        endIdx = startIdx+nMissing-1;
        if endIdx>length(negEgs)
            nMissing = endIdx - length(negEgs);
            endIdx = length(negEgs);
        else
            nMissing = 0;
        end
        negEgsSub = [negEgsSub; negEgs(startIdx:endIdx)];
        startIdx = 1;
    end

    selectedDocs = [posEgs' negEgsSub'];
    selectedDocs = selectedDocs(1:batchSize);

    %save as file to be read by classify and select
    selectedDocs = sort(selectedDocs);
    selectionMap = sparse(1,max(selectedDocs));
    for i=1:length(selectedDocs)
        selectionMap(selectedDocs(i)) = i;
    end
   
    save(selectedDocsFile, 'selectedDocs', 'selectionMap');

    Xfeat = dlmread(featureFile); %should produce Xfeat
    nDocs = length(selectedDocs);    
    
    if size(Xfeat,2)==3 && max(max(Xfeat))>1 %sparse format - rejig this to be a matrix

        nAllDocs = Xfeat(1,1);
        nFeat = Xfeat(2,1);
        Xfeat = Xfeat(3:end,:);

        %add one to indexes since matlab indexes from 1
        Xfeat = sparse(Xfeat(:,1)+1, Xfeat(:,2)+1, Xfeat(:,3), nAllDocs,nFeat); 
    else
        nAllDocs = size(Xfeat,1);
        nFeat = size(Xfeat,2);
    end

    Xfeat = full(Xfeat(selectedDocs,:));

    options = statset('UseParallel', 'always');
    if nRandomBootstrap > 0
        cIdx = kmeans(Xfeat, nRandomBootstrap, 'EmptyAction', 'singleton', 'options', options);
        [cIds, docIdxs] = unique(cIdx);
        bootstrapSet = docIdxs;
    else
        bootstrapSet = [];
    end
    
%     display('Getting the workers to repeat the same tasks initially');
%     workersToRequest = repmat(1:nWorkers, 1, nRandomBootstrap);
%     workersToRequest = workersToRequest(1:(nRandomBootstrap*nWorkers));
%     samplesToRequest = reshape(repmat(selectedDocs(bootstrapSet), nWorkers, 1), 1, nWorkers*nRandomBootstrap);    
    
    workersToRequest = repmat(1:nWorkers, 1, ceil(nRandomBootstrap/nWorkers));
    workersToRequest = workersToRequest(1:nRandomBootstrap);
    samplesToRequest = selectedDocs(bootstrapSet);
    nLabels = 0;

    clear agents;
    clear initialLabels
    initAgents;
    
    [initialLabels, responders] = getSimResponses(agents, workersToRequest, ...
            samplesToRequest, [], qRels(:,[1; chosenIdx]) );
        
    %     [~, initialLabels(:,5)] = find(qRels(initialLabels(:,2)
    save(bootstrapFile, 'bootstrapSet', 'initialLabels');
    
end
clear initialLabels

