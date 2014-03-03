%Simulate the hiring process on static data collected using the old method

%Features for each file
featureFile = '/homes/49/edwin/data/aaai/subset_features.csv';

%Maps original file IDs to the IDs in the feature file
fileMapFile = '/homes/49/edwin/matlab/data/aaai/fileMap_sorted.csv';

%Label file for crowdsourced labels
simulation = false;
%The real labels
sourceLabelFile = '/homes/49/edwin/data/aaai/test_hiring/crowd_class_labels_IJCAI_5.csv';
labelFile = '/homes/49/edwin/data/aaai/test_hiring/simCrowdLabels.csv';

%For writing requests for each round
outputRequestFile = '/homes/49/edwin/data/aaai/test_hiring/outputRequestFile.csv';
excludedWorkersFile = '/homes/49/edwin/data/aaai/test_hiring/excludedWorkers.csv';

%For writing results in the correct format for the test pairs
outputResultFile = '/homes/49/edwin/data/aaai/test_hiring/outputResultFile.csv';

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

nRounds = 400;

currentLabels = zeros(0,5);

% Run the TEST ALGORITHM --------------------------------------------------
nToTry = 1;    
    
crowdTrust = [10 11]
    
excludedWorkers = [];
aucs = zeros(size(qRels,2), nRounds);
lam = zeros(1,nRounds);

workerQueue = 0;

% for f=1:nRounds 
f = 0;
x = [];

sourceLabels = dlmread(sourceLabelFile);
currentLabels = sourceLabels(1:100,:);
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
    
    while ~isempty(workerQueue)
        workerToAllocate = workerQueue(1);
        if length(workerQueue)>1
            workerQueue = workerQueue(2:end);
        else
            workerQueue = [];
        end
        [s, w, excludedWorkers, P, combiner, currentResps, origTsorted] = ...
            classifyAndSelectDoc( nToTry, workerToAllocate, featureFile, labelFile, ...
            excludedWorkersFile, outputRequestFile, excludedWorkersFile, topicDocNoPairsFile, crowdTrust);
    end
         
    binaryMethod = 'roc proportions posterior likelihood';
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
    if f<nRounds && x(end)<length(sourceLabels)-4
        newLabels = sourceLabels(x(end)+1:x(end)+5, :);
        if ~isempty(excludedWorkers)
            newLabels( ismember(newLabels(:,1), excludedWorkers),:) = [];
        end
        workerQueue = [workerQueue newLabels(:,1)'];
        currentLabels = [currentLabels; newLabels];
    end
    fclose('all');
end

% lam(f)
aucs(:,f)'