clear
selectMethod = 'none';    
featureFile = '/homes/49/edwin/data/trec/finalOutput/MatrixTopic2000/Topic-Matrix/matrices_2000_topic.txt';
fileMapFile = '/homes/49/edwin/data/trec/finalOutput/MatrixTopic2000/Topic-Matrix/fileMap.txt';
labelFile = '/homes/49/edwin/data/trec/finalOutput/crowd_class_labels7.csv';
outputRequestFile = '/homes/49/edwin/data/trec/finalOutput/outputRequestFile';
outputResultFile = '/homes/49/edwin/data/trec/posterOutput/output_%s.txt';
topicDocNoPairsFile = '/homes/49/edwin/data/trec/finalOutput/trec2012-crowdsourcing-text-task.topic-docnos.txt';

if ~exist('qrels','var')
    qrels = [];
end
pairsResSet = {};
aucs = {};


runName = 'Orc2Stage';
display(runName)

classifierMethod = 'two-stage';
crowdTrust = [4 16];
confWeight = 1;

[samplesToRequest, P, binaryRes, pairsRes] = classifyAndSelectDoc( 0, 0, ...
    featureFile, labelFile, outputRequestFile, outputResultFile, ...
    topicDocNoPairsFile, fileMapFile, [], [], true, selectMethod, ...
    classifierMethod, confWeight, crowdTrust, runName);
pairsResSet{6} = pairsRes;
[qrels aucs{6}] = evaluateTrecResults(pairsResSet{6}, fileMapFile, qrels, false);
%%%%%%%%%%%%%%

runName = 'Orc2G';
display(runName)

classifierMethod = 'two-stage-G';
crowdTrust = [4 16];
confWeight = 1;

[samplesToRequest, P, binaryRes, pairsRes] = classifyAndSelectDoc( 0, 0, ...
    featureFile, labelFile, outputRequestFile, outputResultFile, ...
    topicDocNoPairsFile, fileMapFile, [], [], true, selectMethod, ...
    classifierMethod, confWeight, crowdTrust, runName);
pairsResSet{1} = pairsRes;
[qrels aucs{1}] = evaluateTrecResults(pairsResSet{1}, fileMapFile, qrels, false);
%%%%%%%%%%%%%

% runName = 'Orc2GUL';
% display(runName)
% 
% classifierMethod = 'two-stage-G-inclLabels';
% crowdTrust = [4 16];
% confWeight = 1;
% 
% [samplesToRequest, P, binaryRes, pairsRes] = classifyAndSelectDoc( 0, 0, ...
%     featureFile, labelFile, outputRequestFile, outputResultFile, ...
%     topicDocNoPairsFile, fileMapFile, [], [], true, selectMethod, ...
%     classifierMethod, confWeight, crowdTrust, runName);
% pairsResSet{2} = pairsRes;
% [qrels aucs{2}] = evaluateTrecResults(pairsResSet{2}, fileMapFile, qrels, false);
% %%%%%%%%%%%%%%
% 
% runName = 'Orc2GULConf';
% display(runName)
% 
% classifierMethod = 'two-stage-G-inclLabels';
% crowdTrust = [4 16];
% confWeight = 'learn';
% 
% [samplesToRequest, P, binaryRes, pairsRes] = classifyAndSelectDoc( 0, 0, ...
%     featureFile, labelFile, outputRequestFile, outputResultFile, ...
%     topicDocNoPairsFile, fileMapFile, [], [], true, selectMethod, ...
%     classifierMethod, confWeight, crowdTrust, runName);
% pairsResSet{3} = pairsRes;
% [qrels aucs{3}] = evaluateTrecResults(pairsResSet{3}, fileMapFile, qrels, false);
% %%%%%%%%%%%%%

% %%%%%And so it begins...
% runName = 'OrcVB1';
% display(runName)
% 
% classifierMethod = 'VB-uncertLabels';
% crowdTrust = [20 80];
% confWeight = 1;
% 
% [samplesToRequest, P, binaryRes, pairsRes] = classifyAndSelectDoc( 0, 0, ...
%     featureFile, labelFile, outputRequestFile, outputResultFile, ...
%     topicDocNoPairsFile, fileMapFile, [], [], true, selectMethod, ...
%     classifierMethod, confWeight, crowdTrust, runName);
% qrels = evaluateTrecResults(pairsRes, fileMapFile, qrels);
% %%%%%%%%%%%%%%%
% 
% runName = 'OrcVB1Conf';
% display(runName)
% 
% classifierMethod = 'VB-uncertLabels';
% crowdTrust = [20 80];
% confWeight = 'learn';
% 
% [samplesToRequest, P, binaryRes, pairsRes] = classifyAndSelectDoc( 0, 0, ...
%     featureFile, labelFile, outputRequestFile, outputResultFile, ...
%     topicDocNoPairsFile, fileMapFile, [], [], true, selectMethod, ...
%     classifierMethod, confWeight, crowdTrust, runName);
% qrels = evaluateTrecResults(pairsRes, fileMapFile, qrels);
% %%%%%%%%%%%%%%%
% 
% runName = 'OrcVBW80'; %PRIMARY SUBMISSION: these labels are the wrong way round in the submission
% display(runName)
% 
% classifierMethod = 'VB-workerUncert';
% crowdTrust = [20 80];
% confWeight = 1;
% 
% [samplesToRequest, P, binaryRes, pairsRes] = classifyAndSelectDoc( 0, 0, ...
%     featureFile, labelFile, outputRequestFile, outputResultFile, ...
%     topicDocNoPairsFile, fileMapFile, [], [], true, selectMethod, ...
%     classifierMethod, confWeight, crowdTrust, runName);
% 
% qrels = evaluateTrecResults(pairsRes, fileMapFile, qrels, false);
% %%%%%%%%%%%%%%%
% 
% runName = 'OrcVBW80Conf'; %these labels are the wrong way round in the submission
% display(runName)
% 
% classifierMethod = 'VB-workerUncert';
% crowdTrust = [20 80];
% confWeight = 'learn';
% 
% [samplesToRequest, P, binaryRes, pairsRes] = classifyAndSelectDoc( 0, 0, ...
%     featureFile, labelFile, outputRequestFile, outputResultFile, ...
%     topicDocNoPairsFile, fileMapFile, [], [], true, selectMethod, ...
%     classifierMethod, confWeight, crowdTrust, runName);
% qrels = evaluateTrecResults(pairsRes, fileMapFile, qrels, false);
% % %%%%%%%%%%%%%%%
% 
% runName = 'OrcVBW9Conf';
% display(runName)
% 
% classifierMethod = 'VB-workerUncert';
% crowdTrust = [1 9];
% confWeight = 'learn';
% 
% [samplesToRequest, P, binaryRes, pairsRes] = classifyAndSelectDoc( 0, 0, ...
%     featureFile, labelFile, outputRequestFile, outputResultFile, ...
%     topicDocNoPairsFile, fileMapFile, [], [], true, selectMethod, ...
%     classifierMethod, confWeight, crowdTrust, runName);
% qrels = evaluateTrecResults(pairsRes, fileMapFile, qrels, false);
% %%%%%%%%%%%%%%

%THIS WAS NOT ONE OF OUR SUBMITTED RUNS AS ONLY 10 WERE ALLOWED
% runName = 'OrcVBW16';
% display(runName)
% 
% classifierMethod = 'VB-workerUncert';
% crowdTrust = [4 16];
% confWeight = 1;
% 
% [samplesToRequest, P, binaryRes, pairsRes] = classifyAndSelectDoc( 0, 0, ...
%     featureFile, labelFile, outputRequestFile, outputResultFile, ...
%     topicDocNoPairsFile, fileMapFile, [], [], true, selectMethod, ...
%     classifierMethod, confWeight, crowdTrust, runName);
% pairsResSet{4} = pairsRes;
% [qrels aucs{4}] = evaluateTrecResults(pairsResSet{4}, fileMapFile, qrels, false);

%%%%%%%%%%%%%

runName = 'OrcVBW16Conf';
display(runName)

classifierMethod = 'VB-workerUncert';
crowdTrust = [4 16];
confWeight = 0.7;

[samplesToRequest, P, binaryRes, pairsRes] = classifyAndSelectDoc( 0, 0, ...
    featureFile, labelFile, outputRequestFile, outputResultFile, ...
    topicDocNoPairsFile, fileMapFile, [], [], true, selectMethod, ...
    classifierMethod, confWeight, crowdTrust, runName);
pairsResSet{5} = pairsRes;
[qrels aucs{5}] = evaluateTrecResults(pairsResSet{5}, fileMapFile, qrels, false);