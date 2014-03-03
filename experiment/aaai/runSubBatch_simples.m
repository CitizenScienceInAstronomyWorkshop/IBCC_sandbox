
clear agents
keepAgents = false;
keepBootstrap = true;

nWorkers = 3;
nDegrade = 2;
nRepeats = 1;%25;
nRandomBootstrap = nWorkers; %one task per worker

batchSize = 2712;%2000;
    
% selectMethod = 'Random';

runTag = 'trec_real_final/';
runRootDir = ['/home/edwin/data/' runTag];
outDir = [runRootDir selectMethod];
if ~exist('outDir', 'dir')
    mkdir(outDir);
end

analDir = ['/home/edwin/results/' runTag];
analysisDataFile = [analDir selectMethod];
if ~exist('analDir', 'dir')
    mkdir(analDir);
end

%     chosenIdx = [3 7 9]';%In order of most examples: 8, 10, 4, 9, 3, 5, 6, 2, 7, 11
%     chosenIdx = [2 5 6]';
chosenIdx = [411 416 417 420 427 432 438 445 446 447]';%[2 3 4 5 6 7 8 9 10 11]';%5;%[2 5 6]';%[8 4 6 10 9 3 5 6 2 7 11]';

if ~exist('aucs_set','var')
    aucs_set = cell(1, nRepeats);
    H_set = cell(1, nRepeats);
    worker_set = cell(1,nRepeats);
end

%This script combines the above. In fact this means we don't need to iterate through the system, just run it with a random bootstrap.
if ~exist('rootDir','var')
    rootDir = '/home/edwin/matlab/data/aaai';
else
    display(['Root dir: ' rootDir]);
end
featureFile = [rootDir '/matrices_2000_topic_norm.txt'];
fileMapFile = [rootDir '/fileMap.txt']; %complete TREC dataset
topicDocNoPairsFile = [rootDir '/trec2012-crowdsourcing-text-task.topic-docnos.txt'];

%Qrels - "ground truth" judgements for testing
qrelFile = [rootDir '/trec-2012-trat-adjudicated-judgments-Oct-11-2012'];
if ~isempty(qrelFile) && ~exist('qRels','var')
    [qRels, qRelsNR] = loadQrels(fileMapFile, qrelFile, topicDocNoPairsFile);
end

qrelFile = [rootDir '/qrels.trec8.adhoc.parts1-5'];
if ~isempty(qrelFile) && ~exist('qRels_old','var')
    [qRels_old, qRelsNR_old] = loadQrels(fileMapFile, qrelFile, topicDocNoPairsFile);
    qRels_old = qRels_old(1:end, 1:size(qRels,2));
%     origTSorted =  [411 416 417 420 427 432 438 445 446 447]';
%     qRels_old = [zeros(size(qRels,1), 1) qRels_old(:, origTSorted)];
end
% subsampleData;

for repeat=1:nRepeats
    display(['repeats: ' num2str(repeat)]);
%     
%     if ~isempty(aucs_set{repeat})
%         display('skipping');
%         continue;
%     end

    runSimDynSub
    
    aucs_set{repeat} = aucs;
    H_set{repeat} = Hs;
    worker_set{repeat} = agents;
       
    save([labelFile '_' num2str(chosenIdx') '_selectMethod' num2str(repeat)], 'currentLabels');
    save(analysisDataFile, 'aucs_set','H_set', 'worker_set');
end

fclose('all');
% 
% chosenIdx = [3 7 9]';
% 
% nRepeats = 5;
% aucs_set2 = cell(1, nRepeats);
% H_set2 = cell(1, nRepeats);
% worker_set2 = cell(1, nRepeats);
% 
% for repeat=1:nRepeats
%     display(['repeats: ' num2str(repeat)]);
% 
%     runSimDynSub
%     
%     aucs_set2{repeat} = aucs;
%     H_set2{repeat} = Hs;
%     worker_set2{repeat} = agents;
%     
%     save([labelFile '_' num2str(chosenIdx') '_selectMethod' num2str(repeat)], 'currentLabels');
%     save(analysisDataFile,  'aucs_set','H_set', 'worker_set', 'aucs_set2','H_set2', 'worker_set2');    
% end
% 
% fclose('all');