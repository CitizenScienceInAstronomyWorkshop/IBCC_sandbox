rootDir = '/homes/49/edwin/data/aaai';

%Features for each file
featureFile = [rootDir '/subset_features.csv'];

%Maps original file IDs to the IDs in the feature file
fileMapFile = [rootDir '/fileMap_sorted.csv'];

%Label file for crowdsourced labels
simulation = true;
%The real labels
labelFile = [rootDir '/sim_dyn/simCrowdLabels.csv'];

%For writing requests for each round
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

nClasses = 11;

%FEATURES -----------------------------------------------------------------
X = dlmread(featureFile, ' ');
nDocs = X(1,1);
nFeatures = X(2,1);

if ~isempty(qrelFile) && ~exist('qRels','var')
    qRels = loadQrels(fileMapFile, qrelFile, topicDocNoPairsFile);
end

nRounds = 800;

if ~exist('nWorkers','var')
    nWorkers = 20;
end

if ~exist('nDegrade','var')
    nDegrade = 0;
end

agents = cell(nWorkers, 1);
for k=1:nWorkers-nDegrade
    agents{k} = SimAgent(false);
end
for k=k:nWorkers
    agents{k} = SimAgent(true, nRounds/(nWorkers));
end
currentLabels = zeros(0,5);

% Run the TEST ALGORITHM --------------------------------------------------
nToTry = 1;    
    
crowdTrust = [50 55]
Nu0 = repmat(1000, 1, nClasses);
%for the features only
Alpha0 = 5 .* ones(nClasses, 2); %0.5 
Alpha0(:,2) = 1; %0.25

excludedWorkers = [];
aucs = zeros(size(qRels,2), nRounds);
aucs_st = zeros(size(qRels,2), nRounds);
lam = zeros(1,nRounds);

workerQueue = 0;

% for f=1:nRounds 
f = 0;
x = [];

%just test a single round to compare the classifier setups not label selection
oneRoundTest = true; 

while size(currentLabels,1)<=nRounds
    f = f+1;
    
    if size(currentLabels,1)>0 && oneRoundTest
        
    dlmwrite(labelFile, currentLabels);
    dlmwrite(excludedWorkersFile, excludedWorkers);
        
    x = [x; size(currentLabels,1)];
    display(['Round ' num2str(f) ', nLabels=' num2str(size(currentLabels,1))]);
    
    if f==nRounds
        writeOutput = false; %true;
    else
        writeOutput = false;
    end

    samplesToRequest = [];
    workersToRequest = [];    
    
    while ~isempty(workerQueue)
        workerToAllocate = workerQueue(1);
        if length(workerQueue)>1
            workerQueue = workerQueue(2:end);
        else
            workerQueue = [];
        end
        [sw, ww, excludedWorkers, P, combiner, currentResps, origTsorted] = ...
            classifyAndSelectDoc( nToTry, workerToAllocate, featureFile, labelFile, ...
            excludedWorkersFile, outputRequestFile, excludedWorkersFile, ...
            topicDocNoPairsFile, crowdTrust, Alpha0, Nu0);
        samplesToRequest = [samplesToRequest sw];
        workersToRequest = [workersToRequest ww];
    end
     
    samplesToRequest(ismember(workersToRequest, excludedWorkers)) = [];
    workersToRequest(ismember(workersToRequest, excludedWorkers)) = [];
    
    %if we've excluded all existing workers, make a new one
    if isempty(workersToRequest)
        agents{nWorkers+1} = SimAgent(true, nRounds/(nWorkers*5));
        nWorkers = nWorkers + 1;
    end    
    
    if length(workersToRequest)<10%f<10 || mod(f,4)==0 || isempty(workersToRequest)
        [s0, w0, excludedWorkers, P, combiner, currentResps, origTsorted] = ...
            classifyAndSelectDoc( nToTry, 0, featureFile, labelFile, ...
                excludedWorkersFile, outputRequestFile, excludedWorkersFile, ...
                topicDocNoPairsFile, crowdTrust, Alpha0, Nu0);   
        samplesToRequest = [samplesToRequest s0];
        workersToRequest = [workersToRequest w0];
%         samplesToRequest = [samplesToRequest sw(length(workerQueue)+1)];
%         workersToRequest = [workersToRequest 0];
    end  
        
    binaryMethod = 'roc proportions posterior likelihood';
%     [binaryRes pairsRes] = writeTestResults( P, combiner, fileMapFile, outputResultFile, ...
%             topicDocNoPairsFile, writeOutput, [], binaryMethod, currentResps, origTsorted );
    
    %EVALUATE LAM AND AUC -------------------------------------------------
    if size(currentLabels,1)>=nRounds
        figure;
        noGraphs = false;
    else
        noGraphs = true;
    end
    aucs_f = evaluateResults(P, [], origTsorted, fileMapFile, ...
        topicDocNoPairsFile, qRels, noGraphs);
    avg_aucs = sum(aucs_f)./(length(aucs_f)-1);
    display(['dyn mean auc: ' num2str(avg_aucs)]);
    aucs(:,f) = aucs_f;
%     lam(f) = lam_f;

    %COMPARE STATIC IBCC WITH CURRENT LABELS ------------------------------
    [~, P_static, ~, ~, origTsorted] = ...
        classifyAndSelectDoc_static( 0, 0, featureFile, labelFile, ...
        outputRequestFile, outputResultFile, topicDocNoPairsFile, fileMapFile, Alpha0, Nu0, ...
        false, 'none', 'VB-workerUncert', 1, crowdTrust);
    aucs_st_f = evaluateResults(P_static, [], origTsorted, fileMapFile, ...
        topicDocNoPairsFile, qRels, noGraphs);
    aucs_st(:,f) = aucs_st_f;
    avg_aucs_st = sum(aucs_st_f)./(length(aucs_st_f)-1);
    display(['sta mean auc: ' num2str(avg_aucs_st)]);
    
    display(['Difference: ' num2str(100*(avg_aucs_st - avg_aucs))]);
    end
    workerQueue = [];    
    
    %Update the set of labels with the requested ones ---------------------
    if f<nRounds && simulation
        %Get the requested labels from the prepared set
        newLabels = [];
        
        display('skipping straight to the end');
        if oneRoundTest
            workersToRequest = randi(nWorkers, 400,1);
            samplesToRequest = randi(nDocs, 400, 1);
        end
        
        for t=1:length(samplesToRequest)
            w = workersToRequest(t);
            i = samplesToRequest(t);
            if w==0
                newWorker = 0;
                while newWorker==0 || ismember(newWorker, excludedWorkers)
                    newWorker = randi(nWorkers,1,1);
                end
                r = agents{newWorker}.getResponse(qRels,i);
                if r~=-1
                    newLabels = [newLabels; newWorker i r 0 0];
%                     workerQueue = [workerQueue newWorker];
                end
            else
                r = agents{w}.getResponse(qRels,i);
                if r~=-1
                    newLabels = [newLabels; w i r 0 0];
%                     workerQueue = [workerQueue w];
                end
            end
            
        end
        currentLabels = [currentLabels; newLabels];
    end
    fclose('all');
end

% lam(f)
aucs(:,f)'