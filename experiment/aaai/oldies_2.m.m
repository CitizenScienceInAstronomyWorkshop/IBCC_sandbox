
%Hypothesis: Intelligent tasking using information gain really only works
%when we have an AUC > 0.5; Otherwise meaningful labels are actually making
%the distributions more uncertain and so aren't valued correctly. Perhaps 
%this is due to IBCC + VB making poor estimates of uncertainty when counts
%are low. Perhaps we can consider biggest changes to labels before summing
%over p(c) rather than taking the total expected information gain, where
%negative gains reduce the utility of a decision but in reality might
%reflect a meaningful change that should have been marked as positive if we
%had a better model/inference algorithm.

clear agents
keepAgents = false;
keepBootstrap = true;

nWorkers = 3;
nDegrade = 2;
nRepeats = 1;
nRandomBootstrap = 30;

batchSize = 500;
    
% selectMethod = 'Random';

runTag = 'hcomp_sf44/';
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

%     chosenIdx = [3 7 9]';%In order of most examples: 8, 10, 4, 9, 3, 5, 6, 2, 7, 11
%     chosenIdx = [2 5 6]';
chosenIdx = [9 10 5]';%6 8 4  10]';

if ~exist('aucs_set','var')
    aucs_set = cell(1, nRepeats);
    H_set = cell(1, nRepeats);
    worker_set = cell(1,nRepeats);
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