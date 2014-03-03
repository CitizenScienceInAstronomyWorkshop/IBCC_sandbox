
clear agents
keepAgents = false;
keepBootstrap = true;

nWorkers = 5;
nDegrade = 0;
nRepeats = 20;
nRandomBootstrap = nWorkers*10; %one task per worker

batchSize = 600;
nToCollect = 600;
% selectMethod = 'Random';

runTag = 'hcomp_woo3/';
runRootDir = ['/homes/49/edwin/data/' runTag];
outDir = [runRootDir selectMethod];
if ~exist('outDir', 'dir')
    mkdir(outDir);
end

analDir = ['/homes/49/edwin/results/' runTag];
analysisDataFile = [analDir selectMethod];
if ~exist('analDir', 'dir')
    mkdir(analDir);
end

alphaFile = [analDir selectMethod '_alpha'];
scoreFile = [analDir selectMethod '_scores'];
unknownFile = [analDir selectMethod '_unknownScores'];

%     chosenIdx = [3 7 9]';%In order of most examples: 8, 10, 4, 9, 3, 5, 6, 2, 7, 11
%     chosenIdx = [2 5 6]';
chosenIdx = 5%[2 5 6]';%[8 4 6 10 9 3 5 6 2 7 11]'; %try 10 3 or 6 (6 was worst in trec results)

if ~exist('aucs_set','var')
    aucs_set = cell(1, nRepeats);
    H_set = cell(1, nRepeats);
    worker_set = cell(1,nRepeats);
    fpr_set = cell(1, nRepeats);
    fnr_set = cell(1, nRepeats);
    brier_set = cell(1, nRepeats);
end

%This script combines the above. In fact this means we don't need to iterate through the system, just run it with a random bootstrap.
if ~exist('rootDir','var')
    rootDir = '/homes/49/edwin/matlab/data/aaai';
else
    display(['Root dir: ' rootDir]);
end
featureFile = [rootDir '/matrices_2000_topic_thresholded.txt'];
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
end
subsampleData;

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
    fpr_set{repeat} = fprs;
    fnr_set{repeat} = fnrs;
    brier_set{repeat} = briers;
    worker_set{repeat} = agents;
       
    save([labelFile '_' num2str(chosenIdx') '_selectMethod' num2str(repeat)], 'currentLabels');
    save(analysisDataFile, 'aucs_set','H_set', 'worker_set', 'fpr_set', 'fnr_set', 'brier_set');
    
    fclose('all');
    close all
end


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