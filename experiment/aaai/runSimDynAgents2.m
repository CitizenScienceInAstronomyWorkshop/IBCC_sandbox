rootDir = '/homes/49/edwin/data/aaai';

%Features for each file
featureFile = [rootDir '/subset_features.csv'];

%Maps original file IDs to the IDs in the feature file
fileMapFile = [rootDir '/fileMap_sorted.csv'];

%Label file for crowdsourced labels
simulation = true;
%The real labels
labelFile = [rootDir '/sim_dyn/simCrowdLabels.csv'];

%For writing requests for each currentRound
outputRequestFile = [rootDir '/sim_dyn/outputRequestFile.csv'];
excludedWorkersFile = [rootDir '/sim_dyn/excludedWorkers.csv'];

%For writing results in the correct format for the test pairs
outputResultFile = [rootDir '/sim_dyn/outputResultFile.csv'];

%Test pairs provided by TREC
% topicDocNoPairsFile = '/homes/49/edwin/data/trec/qRels/topic_docno_pairs_trec7_relevantOnly.csv';
topicDocNoPairsFile = [rootDir '/trec2012-crowdsourcing-text-task.topic-docnos.txt'];

%Qrels - "ground truth" judgements for testing
% qrelFile = '/homes/49/edwin/data/trec/qrels/qrels.trec8.adhoc.parts1-5';
qrelFile = [rootDir '/trec-2012-trat-adjudicated-judgments-Oct-11-2012'];

% delete('centroids.mat');

nClasses = 11;

%FEATURES -----------------------------------------------------------------
X = dlmread(featureFile, ' ');

if ~isempty(qrelFile) && ~exist('qRels','var')
    qRels = loadQrels(fileMapFile, qrelFile, topicDocNoPairsFile);
end

%SUBSAMPLING!
classified = max(qRels,[],2);
dropable = find(classified==0);
dropable = dropable(1:700);
dropable = ismember(X(:,2),dropable);
X(dropable,:) = [];

%have to make sure code ignores blank documents
% -> remove them from X
% -> translate the document IDs from the crowd labels to the desparsified
% IDs. 
% -> Do all this in classify and select
% -> Does nDocs need to change? Where is it used?

nDocs = X(1,1);
nFeatures = X(2,1);

if ~exist('nWorkers','var')
    nWorkers = 2;
end

if ~exist('nDegrade','var')
    nDegrade = 2;
end

nToCollect = 800;
%bootstrap the process by selecting labels at random; 
%can use this if you just want to test the classifiers for one round, 
%not the intelligent tasking
nRandomBootstrap = 400; 

%intelligent tasking sample size. Set to 0 if you want to try all items
nToTry = 0;    

Nu0 = repmat(1000, 1, nClasses);
Nu0(1) = 5000;

%alpha for the workers
crowdTrust = [...
    21 5 5 5 5 5 5 5 5 5 5; ...
    20 6 5 5 5 5 5 5 5 5 5; ...
    20 5 6 5 5 5 5 5 5 5 5; ...
    20 5 5 6 5 5 5 5 5 5 5; ...
    20 5 5 5 6 5 5 5 5 5 5; ...
    20 5 5 5 5 6 5 5 5 5 5; ...
    20 5 5 5 5 5 6 5 5 5 5; ...
    20 5 5 5 5 5 5 6 5 5 5; ...
    20 5 5 5 5 5 5 5 6 5 5; ...
    20 5 5 5 5 5 5 5 5 6 5; ...
    20 5 5 5 5 5 5 5 5 5 6; ...
    ] .* 0.2;

%for the features only: pseudo count for feature not present
Alpha0 = 2 .* ones(nClasses, 2); %0.5 
%pseudo count for feature present
Alpha0(:,2) = 1;%0.25; %0.25

%%%% Initialisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
excludedWorkers = [];
aucs = zeros(size(qRels,2)-1, nToCollect);
aucs_st = zeros(size(qRels,2), nToCollect);
lam = zeros(1,nToCollect);

workerQueue = 0;
nAllocateQueue = 10;

r = 0;
nLabelsEachRound = zeros(1,nToCollect);

noGraphs = true;

workersToRequest = randi(nWorkers, nRandomBootstrap,1);
samplesToRequest = randi(nDocs, nRandomBootstrap, 1);
nLabels = 0;
agents = cell(nWorkers, 1);
for k=1:nWorkers-nDegrade
    agents{k} = SimAgent(k, false);
end
if(numel(k)==0)
    k=1;
end
for k=k:nWorkers
    agents{k} = SimAgent(k, true, 30);
end
currentLabels = zeros(0,5);

%%%%%%%%%% Run! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while nLabels<nToCollect
    r = r+1;
    
    %Update the set of labels with the requested ones ---------------------
    if nLabels<nToCollect && simulation        
        [newLabels, responders] = getSimResponses(agents, workersToRequest, ...
            samplesToRequest, excludedWorkers, qRels);
        currentLabels = [currentLabels; newLabels];
        workerQueue = responders;
        
        [uniqueW, firstIdxs, uniqueIdxs] = unique(workerQueue);
        nAllocateQueue = sparse(1, uniqueIdxs, 1);
        workerQueue = uniqueW;         
    end
    nLabels = size(currentLabels,1);

    dlmwrite(labelFile, currentLabels);
    dlmwrite(excludedWorkersFile, excludedWorkers);

    nLabelsEachRound(r) = nLabels;
    display(['Round ' num2str(r) ', nLabels=' num2str(nLabels)]);

    if nLabels==nToCollect
        writeOutput = false; %true;
    else
        writeOutput = false;
    end
    
    %%% Run Classifier %%%%%%%%%%%%%%% TODO!!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%% allocate workers in the queue %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %select pairs: if we select pairs from a "pool", then one worker will
    %be assigned tasks many times while the others disappear after
    %starvation. So it makes sense to assign to all those with EIG better than
    %unknown workers until we have enough workers, under the assumption
    %that we will need at least n more workers even if the better ones are
    %currently quite good (workers retire before we have enough, workers
    %degrade, multiple workers reduces uncertainty)
    %still leaves us making bespoke selections and automates hire/fire
    %decisions.
    %Alternative is batch of decision requests: say 10 per round. Best
    %worker likely to get lots of tasks,others get replaced - do we limit this?
    %other version is effectively like this but with only one task per
    %person per allocation round. Beware workers would be waiting. So can
    %do the queue as a single batch, make 10 allocations, allowing one
    %worker to have as many allocations as possible. 
    %Problem with original method: (1) we are not selecting between workers 
    %very well - only discarding if they are worse than unknown. (2) not
    %necessarily assigning best person to best task, e.g. if multiple
    %workers have same optimum task. Better to do a small batch for this reason?
    %Original method did batches. Ignore problem of waiting for users,
    %potential for changes in user within a small batch etc., use a greedy
    %method. Ends up focussing on one user? Could that be a good thing
    %provided we have enough time?
    %Process:
    %1. Intialise with 10 unknown workers.
    %2a) try using current queue method to see if searching per worker is
    %possible.
    %2b) Select 20 best pairs; include possibility of new workers.
    %3. <Check this first step selection - is it diverse? If not diverse, 
    %can it waste time on bad workers? Is there a good way to limit
    %tasks per worker? Perhaps we can optimise the selection according to
    %Jensen's inequality, where we reduce risk (minimum information gain?)
    %by diversifying even if this brings down the estimiated Expected IG
    %slightly.
    %>
    %4. Repeat step 2...
    
    %Supply a number of top tasks per iteration. Myopic or greedy selection? 
    %Will it result in only one being chosen? Can have multiple allocated
    %tasks and replace the set each time one is completed. 
    
    nToAllocate = sum(nAllocateQueue);
    samplesToRequest = zeros(1, nToAllocate);
    workersToRequest = zeros(1, nToAllocate);    
    
    while ~isempty(workerQueue)
        %allocate one by one
%         workerToAllocate = workerQueue(1);
%         nToAllocate = nAllocateQueue(1);

        workerToAllocate = workerQueue;
        nK = length(workerToAllocate); %number of workers to allocate at once

        if length(workerQueue)>1
            workerQueue = workerQueue(nK+1:end);
            nAllocateQueue = nAllocateQueue(nK+1:end);
        else
            workerQueue = [];
            nAllocateQueue = [];
        end
        [sw, ww, excludedWorkers, P, combiner, currentResps, origTsorted] = ...
            classifyAndSelectDoc2( nToTry, workerToAllocate, nToAllocate, ...
            featureFile, labelFile, ...
            excludedWorkersFile, outputRequestFile, excludedWorkersFile, ...
            topicDocNoPairsFile, crowdTrust, Alpha0, Nu0);

        newWorkersNeeded = sum(ismember(0,ww));
        
        while newWorkersNeeded > 0
            %replace with new worker
            display('adding a new sim agent');
            agents{nWorkers+1} = SimAgent(true, nToCollect/(nWorkers*5));
            nWorkers = nWorkers + 1;
            newWorkersNeeded = newWorkersNeeded - 1;
        end
%         samplesToRequest = [samplesToRequest sw];
%         workersToRequest = [workersToRequest ww];
        samplesToRequest(1:nToAllocate) = sw;
        workersToRequest(1:nToAllocate) = ww;
    end
    
    %collapse duplicate entries in the queue
    
    %exclude blocked workers
    samplesToRequest(ismember(workersToRequest, excludedWorkers)) = [];
    workersToRequest(ismember(workersToRequest, excludedWorkers)) = [];

    %%%%% Write TREC output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    binaryMethod = 'roc proportions posterior likelihood';
%     [binaryRes pairsRes] = writeTestResults( P, combiner, fileMapFile, outputResultFile, ...
%             topicDocNoPairsFile, writeOutput, [], binaryMethod, currentResps, origTsorted );

    %%%%  EVALUATE LAM AND AUC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nLabels>=nToCollect
        figure;
        noGraphs = false;
    end
    aucs_f = evaluateResults(P, [], origTsorted, fileMapFile, ...
        topicDocNoPairsFile, qRels, noGraphs);
    avg_aucs = sum(aucs_f)./length(aucs_f);
    display(['dyn mean auc: ' num2str(avg_aucs) ': ' num2str(aucs_f')]);
    aucs(:,r) = aucs_f;

%         %COMPARE STATIC IBCC WITH CURRENT LABELS ------------------------------
%         [~, P_static, ~, ~, origTsorted] = ...
%             classifyAndSelectDoc_static( 0, 0, featureFile, labelFile, ...
%             outputRequestFile, outputResultFile, topicDocNoPairsFile, fileMapFile, Alpha0, Nu0, ...
%             false, 'none', 'VB-workerUncert', 1, crowdTrust);
%         aucs_st_f = evaluateResults(P_static, [], origTsorted, fileMapFile, ...
%             topicDocNoPairsFile, qRels, noGraphs);
%         aucs_st(:,r) = aucs_st_f;
%         avg_aucs_st = sum(aucs_st_f)./(length(aucs_st_f)-1);
%         display(['sta mean auc: ' num2str(avg_aucs_st)]);
% 
%         display(['Difference: ' num2str(100*(avg_aucs_st - avg_aucs))]);  

    fclose('all');
end

% lam(r)
aucs(:,r)'