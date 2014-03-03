%Features for each file
featureFile = '/homes/49/edwin/data/aaai/test2/subset_features.csv';

%Maps original file IDs to the IDs in the feature file
fileMapFile = '/homes/49/edwin/matlab/data/aaai/test2/fileMap_sorted.csv';

%Label file for crowdsourced labels
simulation = true;
%The real labels
labelFile = '/homes/49/edwin/data/aaai/test_sim_static/simCrowdLabels.csv';

%For writing requests for each round
outputRequestFile = '/homes/49/edwin/data/aaai/test_sim_static/outputRequestFile.csv';
excludedWorkersFile = '/homes/49/edwin/data/aaai/test_sim_static/excludedWorkers.csv';

%For writing results in the correct format for the test pairs
outputResultFile = '/homes/49/edwin/data/aaai/test_sim_static/outputResultFile.csv';

%Test pairs provided by TREC
% topicDocNoPairsFile = '/homes/49/edwin/data/trec/qRels/topic_docno_pairs_trec7_relevantOnly.csv';
topicDocNoPairsFile = '/homes/49/edwin/data/trec/trec8_sample/trec2012-crowdsourcing-text-task.topic-docnos.txt';

%Qrels - "ground truth" judgements for testing
% qrelFile = '/homes/49/edwin/data/trec/qrels/qrels.trec8.adhoc.parts1-5';
qrelFile = '/homes/49/edwin/data/trec/trec-2012-trat-adjudicated-judgments-Oct-11-2012';

nClasses = 10;

%FEATURES -----------------------------------------------------------------
X = dlmread(featureFile, ' ');
nDocs = X(1,1);
nFeatures = X(2,1);

if ~isempty(qrelFile) && ~exist('qRels','var')
    qRels = loadQrels(fileMapFile, qrelFile, topicDocNoPairsFile);
end

nRounds = 800;

agents = cell(nWorkers, 1);
for k=1:nWorkers-nDegrade
    agents{k} = SimAgent(false);
end
for k=k:nWorkers
    agents{k} = SimAgent(true, nRounds/(nWorkers));
end
currentLabels = [1 1 agents{1}.getResponse(qRels,1) 0 0];

% Run the TEST ALGORITHM --------------------------------------------------
nToTry = 100;
nRequests = 20;

crowdTrust = [500 510]
    
excludedWorkers = [];
aucs = zeros(size(qRels,2), nRounds);
lam = zeros(1,nRounds);

workerQueue = 0;

% for f=1:nRounds 
f = 0;
x = [];
while size(currentLabels,1)<nRounds
    f = f+1;
    dlmwrite(labelFile, currentLabels);
    dlmwrite(excludedWorkersFile, excludedWorkers);
        
    x = [x; size(currentLabels,1)];    
    display(['Round ' num2str(f) ', nLabels=' num2str(size(currentLabels,1))]);
    
    if f==nRounds
        writeOutput = false; %true;
    else
        writeOutput = false;
    end   
        
    [samplesToRequest, P, ~, ~, origTsorted] = ...
            classifyAndSelectDoc( nToTry, nRequests, featureFile, labelFile, ...
            outputRequestFile, outputResultFile, topicDocNoPairsFile, fileMapFile, [], [], ...
            false, 'uncert', 'VB-workerUncert', 1, crowdTrust);
    workersToRequest = randi(nWorkers,1,50);
    
    
%     binaryMethod = 'roc proportions posterior likelihood';
%     [binaryRes pairsRes] = writeTestResults( P, combiner, fileMapFile, outputResultFile, ...
%             topicDocNoPairsFile, writeOutput, [], binaryMethod, currentResps, origTsorted );
    
    %EVALUATE LAM AND AUC -------------------------------------------------
    if f==nRounds
        figure;
        noGraphs = false;
    else
        noGraphs = true;
    end
    [aucs_f] = evaluateResults(P, [], origTsorted, fileMapFile, topicDocNoPairsFile, qRels, noGraphs);
    display(['mean auc: ' num2str(sum(aucs_f)./(length(aucs_f)-1))]);
    aucs(:,f) = aucs_f;
%     lam(f) = lam_f;

    %Update the set of labels with the requested ones ---------------------
    if f<nRounds && simulation
        %Get the requested labels from the prepared set
        newLabels = [];
        
        for t=1:length(samplesToRequest)
            w = workersToRequest(t);
            i = samplesToRequest(t);
            if w==0
                newWorker = randi(nWorkers,1,1);
                r = agents{newWorker}.getResponse(qRels,i);
                if r==-1
                    r = 1;
                end
                newLabels = [newLabels; newWorker i r 0 0];
                workerQueue = [workerQueue newWorker];
                
            else
                r = agents{w}.getResponse(qRels,i);
                if r==-1
                    r = 1;
                end
                newLabels = [newLabels; w i r 0 0];
                workerQueue = [workerQueue w];
                
            end
            
        end
        currentLabels = [currentLabels; newLabels];
    end
    fclose('all');
end

% lam(f)
aucs(:,f)'