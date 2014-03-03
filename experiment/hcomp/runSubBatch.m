
clear agents
keepAgents = false;
keepBootstrap = true;

nRepeats = 10;
nRandomBootstrap = 10;

batchSize = 500;
    
bootstrapSet = {};
for r=1:nRepeats
    bootstrapSet{r} = (1:nRandomBootstrap) + ((r-1)*nRandomBootstrap);
end

% selectMethod = 'Random';
outDir = ['/homes/49/edwin/data/hcomp_sf4/' selectMethod];
if ~exist('outDir', 'dir')
    mkdir(outDir);
end

analDir = '/homes/49/edwin/results/hcomp_sf4/';
analysisDataFile = [analDir selectMethod];
if ~exist('analDir', 'dir')
    mkdir(analDir);
end

if ~exist('chosenIdx','var')
%     chosenIdx = [3 7 9]';

    chosenIdx = [2 5 6]';
end

aucs_set = cell(1, nRepeats);
H_set = cell(1, nRepeats);
worker_set = cell(1,nRepeats);

for repeat=1:nRepeats
    display(['repeats: ' num2str(repeat)]);

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