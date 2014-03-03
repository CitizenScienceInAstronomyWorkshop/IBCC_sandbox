%Features for each file
%  featureFile = '/homes/49/edwin/data/trec/finalOutput/MatrixTopic2000/Topic-Matrix/matrices_2000_topic.txt';
featureFile = '/homes/49/edwin/data/aaai/test2/subset_features.csv';
% featureFile = '/homes/49/edwin/matlab/data/aaai/subset/subset_features_2.csv';

%Maps original file IDs to the IDs in the feature file
%  fileMapFile = '/homes/49/edwin/data/trec/finalOutput/MatrixTopic2000/Topic-Matrix/fileMap.txt';
%  fileMapFile = '/homes/49/edwin/matlab/data/aaai/subset/fileMap_2.txt';
fileMapFile = '/homes/49/edwin/matlab/data/aaai/test2/fileMap_sorted.csv';

%Label file for crowdsourced labels
% labelFile = '/homes/49/edwin/data/aaai/testOutput/labelFile.csv'; 
% simulation = true; %set to false if using real labels
%The real labels
%  labelFile = '/homes/49/edwin/data/trec/finalOutput/crowd_class_labels7.csv';
labelFile = '/homes/49/edwin/data/aaai/test2/crowd_labels_5.csv';

simulation = false;

%For writing requests for each round
outputRequestFile = '/homes/49/edwin/data/aaai/testOutput/outputRequestFile.csv';
excludedWorkersFile = '/homes/49/edwin/data/aaai/testOutput/excludedWorkers.csv';

%For writing results in the correct format for the test pairs
outputResultFile = '/homes/49/edwin/data/aaai/testOutput/outputResultFile.csv';

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

if simulation && ~exist('crowdLabels','var')
    nInitialLabels = 1000;
    nWorkers = 50;
    [crowdLabels, currentLabels, qRels] = genCrowdLabels(fileMapFile, ...
        topicDocNoPairsFile, qRels, nInitialLabels, nWorkers, nClasses);
end

% Run the TEST ALGORITHM --------------------------------------------------
if simulation
    nRounds = 10;
    workerToAllocate = 0;
    nToTry = 10;    
else
    nRounds = 1;
    workerToAllocate = 10000;
    nToTry = 1;
end
    
crowdTrust = [500 510]
    
excludedWorkers = [];
aucs = zeros(size(qRels,2), nRounds);
lam = zeros(1,nRounds);

for f=1:nRounds 
    if simulation
        dlmwrite(labelFile, currentLabels);
        dlmwrite(excludedWorkersFile, excludedWorkers);
    end
        
    display(num2str(f));
    
    if f==nRounds
        writeOutput = false; %true;
    else
        writeOutput = false;
    end
  
    [samplesToRequest, workersToRequest, excludedWorkers, P, combiner, currentResps, origTsorted] = ...
        classifyAndSelectDoc( nToTry, workerToAllocate, featureFile, labelFile, ...
        excludedWorkersFile, outputRequestFile, excludedWorkersFile, topicDocNoPairsFile, crowdTrust);

    binaryMethod = 'roc proportions posterior likelihood';
    [binaryRes pairsRes] = writeTestResults( P, combiner, fileMapFile, outputResultFile, ...
            topicDocNoPairsFile, writeOutput, [], binaryMethod, currentResps, origTsorted );
    
    %EVALUATE LAM AND AUC -------------------------------------------------
    if f==nRounds
        figure;
        noGraphs = false;
    else
        noGraphs = true;
    end
    [aucs_f, lam_f] = evaluateResults(P, binaryRes, origTsorted, fileMapFile, topicDocNoPairsFile, qRels, noGraphs);
    display(['mean auc: ' num2str(sum(aucs_f)./(length(aucs_f)-1))]);
    aucs(:,f) = aucs_f;
    lam(f) = lam_f;

    %Update the set of labels with the requested ones ---------------------
    if f<nRounds && simulation
        %Get the requested labels from the prepared set
        tmp = crowdLabels(ismember(crowdLabels(:,2),samplesToRequest), :);
        tmp2 = [];
        for t=1:length(samplesToRequest)
            for w=workersToRequest{t}
                if w==0
                    tmp2 = [tmp2; tmp(tmp(:,2)==samplesToRequest(t))];
                else
                    tmp2 = [tmp2; tmp(tmp(:,1)==w,:)];
                end
            end
        end
        currentLabels = [currentLabels; tmp2];
    end
end

lam(f)
aucs(:,f)'