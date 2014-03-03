%No Hiring and Firing - keep using poor workers; can still assign them the most suitable task, don't compare with the unknown worker; start this process by spawning a number of tasks open to all, add more in if stagnation (timeouts on tasks?)
%No Active Selection - give tasks to workers at random.

if ~exist('repeat','var')
    repeat = 1;
end

%This script combines the above. In fact this means we don't need to iterate through the system, just run it with a random bootstrap.
if ~exist('rootDir','var')
    rootDir = '/homes/49/edwin/matlab/data/aaai';
else
    display(['Root dir: ' rootDir]);
end

if ~exist('outDir','var')
    outDir = [rootDir '/simDyn'];
else
    display(['Output dir: ' outDir]);
end

if ~exist('runRootDir','var')
    runRootDir = outDir;
else
    display(['Run-specific root dir: ' rootDir]);
end

%Features for each file
% featureFile = [rootDir '/subset_features.csv']; %IJCAI subset
featureFile = [rootDir '/matrices_2000_topic_thresholded.txt']; %complete TREC dataset. Thresholded features with value < 0.01 for speed (around half all entries).
% featureFile = [rootDir '/synth_feat_3_7_10.txt']; %complete TREC dataset
% if exist('chosenIdx','var') && sum(chosenIdx==[2 5 6]')==3
% %     featureFile = [rootDir '/trim_feat.txt']; % for 2 5 6
%     featureFile = [rootDir '/trim_feat_200docs.txt']; % for 2 5 6
% elseif exist('chosenIdx','var') && sum(chosenIdx==[3 7 9]')==3
%     featureFile = [rootDir '/trim_feat379.txt']; % for 3 7 9
% end


%Maps original file IDs to the IDs in the feature file
% fileMapFile = [rootDir '/fileMap_sorted.csv']; %IJCAI subset
fileMapFile = [rootDir '/fileMap.txt']; %complete TREC dataset

%Label file for crowdsourced labels
simulation = true;
%The real labels
labelFile = [outDir '/simCrowdLabels.csv'];

%For writing requests for each currentRound
outputRequestFile = [outDir '/outputRequestFile.csv'];
excludedWorkersFile = [outDir '/excludedWorkers.csv'];

%For writing results in the correct format for the test pairs
outputResultFile = [outDir '/outputResultFile.csv'];

%Test pairs provided by TREC
% topicDocNoPairsFile = '/homes/49/edwin/data/trec/qRels/topic_docno_pairs_trec7_relevantOnly.csv';
topicDocNoPairsFile = [rootDir '/trec2012-crowdsourcing-text-task.topic-docnos.txt'];

%Qrels - "ground truth" judgements for testing
qrelFile = [rootDir '/trec-2012-trat-adjudicated-judgments-Oct-11-2012'];

%File listing the subsample of documents we are actually testing on
selectedDocsFile = [runRootDir '/selectedDocs_' num2str(repeat) '.mat'];
bootstrapFile = [runRootDir '/bootstrap_CrowdLabels' num2str(repeat) '.mat'];

load(selectedDocsFile);
load(bootstrapFile);

centroidsFile = [outDir '/centroids.mat'];
if exist(centroidsFile,'file')
    delete(centroidsFile);
end

display('use a conf matrix because performance over the very common negative class invariably looks good, but does not say much about ability to pick out the rare positive class');


%FEATURES -----------------------------------------------------------------
X = dlmread(featureFile);

if ~isempty(qrelFile) && ~exist('qRels','var')
    [qRels, qRelsNR] = loadQrels(fileMapFile, qrelFile, topicDocNoPairsFile);
end

qrelFile = [rootDir '/qrels.trec8.adhoc.parts1-5'];
if ~isempty(qrelFile) && ~exist('qRels_old','var')
    [qRels_old, qRelsNR_old] = loadQrels(fileMapFile, qrelFile, topicDocNoPairsFile);
end

% SUBSAMPLING -------------------------------------------------------------

nClasses = length(chosenIdx)+1;

%have to make sure code ignores blank documents
% -> remove them from X
% -> translate the document IDs from the crowd labels to the desparsified
% IDs.
% -> Do all this in classify and select
% -> Does nDocs need to change? Where is it used?

%BASIC PARAMS -------------------------------------------------------------

nDocs = length(selectedDocs);%X(1,1);
nFeatures = X(2,1);
clear X

%Is this enough to detect boredom?
nToCollect = 501; % 400
%bootstrap the process by selecting labels at random;
%can use this if you just want to test the classifiers for one round,
%not the intelligent tasking
if ~exist('keepBootstrap','var') || ~exist('nRandomBootstrap','var') || keepBootstrap==false
    nRandomBootstrap = 5;
end

if ~exist('selectMethod','var')
    selectMethod = 'HFAS';
end

%intelligent tasking sample size. Set to 0 if you want to try all items
nToTry = 0;

Nu0 = repmat(100, 1, nClasses);
Nu0(1) = 100;

%alpha for the workers
% crowdTrust = [...
%     21 5 5 5 5 5 5 5 5 5 5; ...
%     20 6 5 5 5 5 5 5 5 5 5; ...
%     20 5 6 5 5 5 5 5 5 5 5; ...
%     20 5 5 6 5 5 5 5 5 5 5; ...
%     20 5 5 5 6 5 5 5 5 5 5; ...
%     20 5 5 5 5 6 5 5 5 5 5; ...
%     20 5 5 5 5 5 6 5 5 5 5; ...
%     20 5 5 5 5 5 5 6 5 5 5; ...
%     20 5 5 5 5 5 5 5 6 5 5; ...
%     20 5 5 5 5 5 5 5 5 6 5; ...
%     20 5 5 5 5 5 5 5 5 5 6; ...
%     ] .* 0.2;

crowdTrust = ones(nClasses);
crowdTrust( sub2ind(size(crowdTrust), 1:nClasses, 1:nClasses) ) = 3;
% crowdTrust(2:end,1) = 3;
% crowdTrust(1,2:end) = (3+nClasses-2)/(nClasses-1);
crowdTrust = (crowdTrust) + 9;
% crowdTrust = crowdTrust ./ 2;%max(max(crowdTrust));
% crowdTrust = crowdTrust + 1;

%for the features only: pseudo count for feature not present
Alpha0 = 1 .* ones(nClasses, 2); %0.5
%pseudo count for feature present
Alpha0(:,2) = 10;%0.5;%990;%0.25;%0.25; %0.25

%%%% Initialisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
excludedWorkers = [];
aucs = zeros(nClasses, nToCollect);
aucs_st = zeros(nClasses, nToCollect);
Hs = zeros(1, nToCollect); %the entropy of the target labels
lam = zeros(1,nToCollect);

workerQueue = 0;
nAllocateQueue = 1;
allocationCap = 1; %don't allocate more than this number to any one worker

r = 0;
prevNLabels = 0;
nLabelsEachRound = zeros(1,nToCollect);

noGraphs = true;

initAgents;

currentLabels = initialLabels;
nLabels = size(currentLabels,1);
    
dlmwrite(labelFile, currentLabels);
[~, ~, ~, P, combiner] = classifyAndSelectSub( selectMethod, 0, 0, 0, ...
            featureFile, labelFile, selectedDocsFile,...
            excludedWorkersFile, centroidsFile, outputRequestFile, excludedWorkersFile,...
            chosenIdx, crowdTrust, Alpha0, Nu0);
aucs_f = zeros(size(P,2),1);
for j=2:size(P,2)
    testLabels = qRels(selectedDocs,chosenIdx(j-1))';
    testIdxs = selectedDocs;%(testLabels~=0);
%         testLabels = testLabels(testLabels~=0);
    testLabels(testLabels==-1) = 0;
    auc = graphs.ClassifierPerformanceGraph.drawRoc(P(testIdxs,j)', ...
        testLabels, {'topic'}, false, noGraphs, false);
    aucs_f(j) = auc;
end

avg_aucs = sum(aucs_f)./(length(aucs_f)-1);
display(['before intelligent tasking: mean auc: ' num2str(avg_aucs) ': ' num2str(aucs_f')]);        
        

%%%%%%%%%% Run! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while nLabels<nToCollect
    r = r+1;

    %Update the set of labels with the requested ones ---------------------
    if nLabels<nToCollect && simulation && (r>1 || ~exist('initialLabels','var'))
        [newLabels, responders] = getSimResponses(agents, workersToRequest, ...
            samplesToRequest, excludedWorkers, qRels(:,[1; chosenIdx]) );

%         newLabels(:,3) = chosenMap(newLabels(:,3));
%         newLabels(newLabels(:,3)==1,3) = 0;

        currentLabels = [currentLabels; newLabels];
%         if r==1
%             initialLabels = currentLabels;
%         end
    end

    if simulation && nLabels<nToCollect
        if r==1
            responders = currentLabels(:,1);
        end

        workerQueue = responders;

        [uniqueW, firstIdxs, uniqueIdxs] = unique(workerQueue);
        nAllocateQueue = sparse(1, uniqueIdxs, 1);
        nAllocateQueue(nAllocateQueue>allocationCap) = allocationCap;
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
    lastIdx = 1;
    while nLabels<nToCollect && ~isempty(workerQueue)
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
        [sw, ww, excludedWorkers, P, combiner, currentResps, origTsorted, logJoint] = ...
            classifyAndSelectSub( selectMethod, nToTry, workerToAllocate, nToAllocate, ...
            featureFile, labelFile, selectedDocsFile,...
            excludedWorkersFile, centroidsFile, outputRequestFile, excludedWorkersFile,...
            chosenIdx, crowdTrust, Alpha0, Nu0);

        newWorkersNeeded = sum(ww==0);

        while newWorkersNeeded > 0
            %replace with new worker
            display(['adding a new sim agent, id=' num2str(nWorkers+1)]);

            workerType = mod(nWorkers+1, 3) + 1;
            skills = [0.9 0.4 0.65]; %hcomp sf4 has 0.4 moved to the end
            skill = skills(workerType);

%             agents{nWorkers+1} = SimAgent(nWorkers+1, true, 10, skill, nClasses, 5);
            agents{nWorkers+1} = SimAgent(nWorkers+1, true, 1, skill, nClasses, 10);
%             agents{nWorkers+1} = SimAgent(nWorkers+1, false, 0, 1/nClasses, nClasses);

            nWorkers = nWorkers + 1;
            newWorkersNeeded = newWorkersNeeded - 1;
            
            newReqIdx = find(ww==0, 1, 'first');
            ww(newReqIdx) = nWorkers;
        end

        for exW=excludedWorkers'
            agents{exW}.fired = true;
        end
%         samplesToRequest = [samplesToRequest sw];
%         workersToRequest = [workersToRequest ww];
        samplesToRequest(lastIdx:nToAllocate) = sw;
        workersToRequest(lastIdx:nToAllocate) = ww;
        lastIdx = nToAllocate+1;
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

        [sw, ww, excludedWorkers, P, combiner, currentResps, origTsorted, logJoint] = ...
        classifyAndSelectSub( selectMethod, 0, 0, 0, ...
        featureFile, labelFile, selectedDocsFile,...
        excludedWorkersFile, centroidsFile, outputRequestFile, excludedWorkersFile, ...
        chosenIdx, crowdTrust, Alpha0, Nu0);
    end
    aucs_f = zeros(nClasses,1);

%     negLabels = (sum(qRelsNR(selectedDocs,chosenIdx),2)==-length(chosenIdx))';
%     auc = graphs.ClassifierPerformanceGraph.drawRoc(P(selectedDocs,1)', ...
%         negLabels, ...
%         {'topic'}, false, noGraphs, false)
%     aucs_f(1) = auc;

    for j=2:size(P,2)
        testLabels = qRels(selectedDocs,chosenIdx(j-1))';
        testIdxs = selectedDocs;%(testLabels~=0);
%         testLabels = testLabels(testLabels~=0);
        testLabels(testLabels==-1) = 0;
        auc = graphs.ClassifierPerformanceGraph.drawRoc(P(testIdxs,j)', ...
            testLabels, {'topic'}, false, noGraphs, false);
        aucs_f(j) = auc;
    end

    avg_aucs = sum(aucs_f)./(length(aucs_f)-1);
    display(['dyn mean auc: ' num2str(avg_aucs) ': ' num2str(aucs_f')]);
    aucs(:,prevNLabels+1:nLabels) = repmat(aucs_f, 1, nLabels-prevNLabels);

    logJoint(isinf(logJoint)) = 0;
    Hs(prevNLabels+1:nLabels) = -sum(sum(combiner.combinedPost .* log(combiner.combinedPost)));
%     Hs(prevNLabels+1:nLabels)  = -sum(sum(exp(logJoint) .* logJoint, 1), 2) + sum(sum(exp(logJoint),1) .* log(sum(exp(logJoint),1)), 2);
%     Hs(prevNLabels+1:nLabels) = -sum(sum(combiner.combinedPost .* logJoint, 1) + log(sum(exp(logJoint),1)), 2);

    display(['H labels: ' num2str(Hs(nLabels)) ]);

    prevNLabels = nLabels;

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
if r>0
    % lam(r)
    aucs(:,r)'
end