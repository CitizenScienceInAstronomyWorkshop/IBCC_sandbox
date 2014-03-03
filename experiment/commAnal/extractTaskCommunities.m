function [gTasks1, Ptasks1, gTasks2, Ptasks2, AlphaTaskComms, PiTaskComms, observedAgents] = ...
    extractTaskCommunities(Alpha, testData, snBaseOutputs, labels, slice, filterMap, sliceStart)
%Alpha - the hyperparameters for the confusion matrices
%testData - data used during testing
%snBaseOutputs - complete dataset; we include tasks from this dataset 

%Find samples that were not included in the test set but were seen by the
%test agents
dataExpd = testData;


size(unique(dataExpd{1}))

if nargin < 7
    sliceStart = 1;
end

if nargin > 4
    dataExpd{1} = dataExpd{1}(sliceStart:slice);
    dataExpd{2} = dataExpd{2}(sliceStart:slice);
    dataExpd{3} = dataExpd{3}(sliceStart:slice);
    
    lastAgent = find(snBaseOutputs{1}==dataExpd{1}(end));
    lastTask = lastAgent(find(snBaseOutputs{2}(lastAgent)==dataExpd{2}(end)));
    lastScore = lastTask(find(snBaseOutputs{3}(lastTask)==dataExpd{3}(end)));
    completeSetIdx = 1:lastScore;
    snBaseOutputs{1} = snBaseOutputs{1}(completeSetIdx);
    snBaseOutputs{2} = snBaseOutputs{2}(completeSetIdx);
    snBaseOutputs{3} = snBaseOutputs{3}(completeSetIdx);
    
    firstAgent = find(snBaseOutputs{1}==dataExpd{1}(1));
    firstTask = firstAgent(find(snBaseOutputs{2}(firstAgent)==dataExpd{2}(1)));
    firstScore = firstTask(find(snBaseOutputs{3}(firstTask)==dataExpd{3}(1)));
    completeSetIdx = firstScore;
    snBaseOutputs{1} = snBaseOutputs{1}(completeSetIdx:end);
    snBaseOutputs{2} = snBaseOutputs{2}(completeSetIdx:end);
    snBaseOutputs{3} = snBaseOutputs{3}(completeSetIdx:end);
end

size(unique(dataExpd{1}))

display(['slice size: ' num2str(numel(dataExpd{1}))]);

extraIdxs = ismember(snBaseOutputs{1}, dataExpd{1}) & ~ismember(snBaseOutputs{2}, dataExpd{2});

dataExpd{1} = [dataExpd{1}; snBaseOutputs{1}(extraIdxs)];
dataExpd{2} = [dataExpd{2}; snBaseOutputs{2}(extraIdxs)];
dataExpd{3} = [dataExpd{3}; snBaseOutputs{3}(extraIdxs)];

size(unique(dataExpd{1}))

display(['slice size: ' num2str(numel(dataExpd{1}))]);

% [N, X] = hist(dataExpd{2}, unique(dataExpd{2}));
% samplesToKeep = N>1;
% samplesToKeep = X(samplesToKeep);
% sampleIdxs = ismember(dataExpd{2},samplesToKeep);
% dataExpd{1} = dataExpd{1}(sampleIdxs);
% dataExpd{2} = dataExpd{2}(sampleIdxs);
% dataExpd{3} = dataExpd{3}(sampleIdxs);


% [N, X] = hist(dataExpd{1}(idsJ1), unique(dataExpd{1}));
% [N2, X2] = hist(dataExpd{1}(idsJ2), unique(dataExpd{1}));
% samplesToKeep1 = N>5; %require at least N classifications
% samplesToKeep2 = N2>5; %require at least N classifications
% samplesToKeep = X(samplesToKeep1 & samplesToKeep2);

[N, X] = hist(dataExpd{1}, unique(dataExpd{1}));
% [N2, X2] = hist(dataExpd{1}(idsJ2), unique(dataExpd{1}));
samplesToKeep = N>10; %require at least N classifications
samplesToKeep = X(samplesToKeep);

sampleIdxs = ismember(dataExpd{1},samplesToKeep);
dataExpd{1} = dataExpd{1}(sampleIdxs);
dataExpd{2} = dataExpd{2}(sampleIdxs);
dataExpd{3} = dataExpd{3}(sampleIdxs);

% idsJ1 = ismember(dataExpd{2}, find(labels==1));
% idsJ2 = ismember(dataExpd{2}, find(labels==2));
% 
% dataExpd{1} = dataExpd{1}(idsJ2);
% dataExpd{2} = dataExpd{2}(idsJ2);
% dataExpd{3} = dataExpd{3}(idsJ2);

size(unique(dataExpd{1}))

data = dataExpd; %snBaseOutputs; %dataExpd; % testData;

if isempty(Alpha) && ~isempty(agentRatings)
    Alpha = agentRatings{1}{3};
end

observedTasks = unique(data{2});
nTasks = numel(observedTasks);

[observedAgents, agentIdxToID, agentIDToIdx] = unique(data{1});
nAgents = length(observedAgents);
if nargin > 5
    AlphaTaskComms = Alpha(:,:,filterMap(observedAgents));
    normTerm = sum(AlphaTaskComms, 2);
    normTerm = repmat(normTerm, 1, size(AlphaTaskComms,2));
    PiTaskComms = AlphaTaskComms ./ normTerm;
else
    AlphaTaskComms = Alpha(:,:,observedAgents);
    normTerm = sum(AlphaTaskComms, 2);
    normTerm = repmat(normTerm, 1, size(AlphaTaskComms,2));
    PiTaskComms = AlphaTaskComms ./ normTerm;    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Adjacency matrix based on which data points were observed
display(['building network for ' num2str(nAgents) ' agents.']);
Adj = zeros(nAgents, nAgents);

totalsByAgent = zeros(nAgents, 1);

for n=observedTasks'
    
    hits = data{2}==n;
    hitsIdxs = agentIDToIdx(hits);
    Adj(hitsIdxs, hitsIdxs) = Adj(hitsIdxs, hitsIdxs) + 1;
    totalsByAgent(hitsIdxs) = totalsByAgent(hitsIdxs) + 1;
end

%normalise
totalsByAgent = repmat(totalsByAgent, 1, nAgents);
Adj = Adj ./ (totalsByAgent + totalsByAgent');

%P: soft membership score table. Rows->Agents. Columns->Communities.
%g: cell array containing lists of agents in each community
display('starting comm detection');
[Ptasks1, gTasks1] = commDetNMF(Adj, 60);
mod = get_modularity(gTasks1, Adj, 'groups');
display(num2str(mod));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Adjacency based on observation and score
% display('building task-score network');
% Adj2 = zeros(nAgents, nAgents);
% 
% scores = [-1 1 3];
% 
% totalsByAgent = zeros(nAgents, 1);
% 
% for n=observedTasks'
%     
%     hits = data{2}==n;
%     for s=scores
%         scoreHits = (data{3}==s) & hits;        
%         hitsIdxs = agentIDToIdx(scoreHits);   
%         Adj2(hitsIdxs, hitsIdxs) = Adj2(hitsIdxs, hitsIdxs) + 1;
%         totalsByAgent(hitsIdxs) = totalsByAgent(hitsIdxs) + 1;        
%     end
% end
% 
% %normalise
% totalsByAgent = repmat(totalsByAgent, 1, nAgents);
% Adj2 = Adj2 ./ (totalsByAgent + totalsByAgent');
% 
% display('starting comm detection');
% [Ptasks2, gTasks2] = commDetNMF(Adj2, 60);
% mod = get_modularity(gTasks2, Adj2, 'groups');
% display(num2str(mod));
Ptasks2 = [];
gTasks2 = [];

end