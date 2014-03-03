%Community detection over the samples/assets

%Find samples that were not included in the test set but were seen by the
%test agents
dataExpd = testData;

extraIdxs = ismember(snBaseOutputs{1}, testData{1}) & ~ismember(snBaseOutputs{2}, testData{2});

dataExpd{1} = [dataExpd{1}; snBaseOutputs{1}(extraIdxs)];
dataExpd{2} = [dataExpd{2}; snBaseOutputs{2}(extraIdxs)];
dataExpd{3} = [dataExpd{3}; snBaseOutputs{3}(extraIdxs)];


[N, X] = hist(dataExpd{2}, unique(dataExpd{2}));
samplesToKeep = N>1;
samplesToKeep = X(samplesToKeep);
sampleIdxs = ismember(dataExpd{2},samplesToKeep);
dataExpd{1} = dataExpd{1}(sampleIdxs);
dataExpd{2} = dataExpd{2}(sampleIdxs);
dataExpd{3} = dataExpd{3}(sampleIdxs);

[N, X] = hist(dataExpd{1}, unique(dataExpd{1}));
samplesToKeep = N>10; %require at least N classifications
samplesToKeep = X(samplesToKeep);
sampleIdxs = ismember(dataExpd{1},samplesToKeep);
dataExpd{1} = dataExpd{1}(sampleIdxs);
dataExpd{2} = dataExpd{2}(sampleIdxs);
dataExpd{3} = dataExpd{3}(sampleIdxs);

data = dataExpd; %snBaseOutputs; %dataExpd; % testData;

[observedTasks sampleIdxToID sampleIDToIdx] = unique(data{2});
nTasks = numel(observedTasks);

[observedAgents, agentIdxToID, agentIDToIdx] = unique(data{1});
nAgents = length(observedAgents);

%adjacency based on common agents
taskAgentMap = sparse(sampleIDToIdx, agentIDToIdx, ones(numel(data{1}),1));

taskAgentMap = taskAgentMap~=0; %flatten values to zeros and ones

% % taskAgentMap = sparse(nTasks, nAgents);
% % mapIdxs = sub2ind(size(taskAgentMap), sampleIDToIdx, agentIDToIdx);
% % taskAgentMap(mapIdxs) = 1; %use this version if we don't count duplicates

Adj = sparse(nTasks, nTasks); %taskAgentMap * taskAgentMap';

for n=1:nTasks
    Adj(n, :) = taskAgentMap * taskAgentMap(n,:)';
end

%Symmetric co-occurrence: normalise by number of samples seen by both
%agents in a pair
%For entries in Adj which are zero we can set totalsBySample to zero and 
%sparsify to save memory
totalsBySample = repmat(sum(taskAgentMap,2), 1, nTasks);
Adj = Adj ./ (totalsBySample + totalsBySample');
[Ps1, gS1] = commDetNMF(Adj, 10);
mod = get_modularity(gS1, Adj, 'groups');
display(num2str(mod));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%adjacency based on common agent-score combination
% taskAgentMap = sparse(sampleIDToIdx, agentIDToIdx, data{3);

taskAgentMap = sparse(nTasks, nAgents);
mapIdxs = sub2ind(size(taskAgentMap), sampleIDToIdx, agentIDToIdx);
taskAgentMap(mapIdxs) = data{3}; %not sure what happens to multiple scores given by same agent

Adj = sparse(nTasks, nTasks);

totalsBySample = zeros(nTasks, 1);

for s=[-1 1 3]
    T = taskAgentMap==s;
    
    scoreMatches = sparse(nTasks, nTasks);
    for n=1:nTasks
        scoreMatches(n, :) = T* T(n,:)';
    end
    
    Adj = Adj + scoreMatches;
    totalsBySample = totalsBySample + sum(T, 2);
end

totalsBySample = repmat(totalsBySample, 1, nTasks);
Adj = Adj ./ (totalsBySample + totalsBySample');
[Ps2, gS2] = commDetNMF(Adj, 10);
mod = get_modularity(gS2, Adj, 'groups');
display(num2str(mod));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%adjacency based on common pi-community members
% Adj = sparse(nTasks, nTasks);
% totalsBySample = zeros(nTasks, 1);
% 
% %Cooccurrence when samples receive the same scores
% for s=[-1 1 3]
%     [p, dataPiComm] = max(P(agentIDToIdx, :), [], 2);
%     T = sparse(sampleIDToIdx, dataPiComm, double(data{3}==s) );
%     
%     piMatches = sparse(nTasks, nTasks);
%     for n=1:nTasks
%         piMatches(n, :) = T* T(n,:)';
%     end    
%     
%     Adj = Adj + piMatches;
%     totalsBySample = totalsBySample + sum(T, 2);
% end
% % [p, dataPiComm] = max(P(agentIDToIdx, :), [], 2);
% % T = sparse(sampleIDToIdx, dataPiComm, ones(numel(data{1}, 1)));
% % Adj = Adj + (T*T');
% 
% totalsBySample = repmat(totalsBySample, 1, nTasks);
% Adj = Adj ./ (totalsBySample + totalsBySample');
% 
% [Ps3, gS3] = commDetNMF(Adj, 10);
% mod = get_modularity(gS3, Adj, 'groups');
% display(num2str(mod));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adj = sparse(nTasks, nTasks);
% 
% %probability of each community given this sample
% [p, dataPiComm] = max(P(agentIDToIdx, :), [], 2);
% gOverlay = sparse(nTasks, numel(g));
% for c=1:numel(g)
%     T = sparse(sampleIDToIdx, 1, P(agentIDToIdx, c));
%     gOverlay(:, c) = T;
%     
% end
% % gOverlay = sparse(sampleIDToIdx, dataPiComm, ones(numel(data{2},1)));
% 
% gOverlayNorm = gOverlay ./ repmat(sum(gOverlay, 2), 1, numel(g));
% 
% %hellinger distance between these
% squaredHD = 1;
% for c=1:numel(g)
%     squaredHD = squaredHD - (gOverlayNorm(:, c)*gOverlayNorm(:,c)').^0.5;
% end
% 
% squaredHD(squaredHD<0) = 0;
% 
% Adj = exp(-squaredHD);
% Adj = 1-squaredHD.^0.5;
% [Ps4, gS4] = commDetNMF(Adj, 50);
% mod = get_modularity(gS4, Adj, 'groups');
% display(num2str(mod));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%try to find communities related to correctness of answers
Adj = sparse(nTasks, nTasks);

% negTasks = labels(observedTasks)==1;
% posTasks = labels(observedTasks)==2;
% 
% for s=[-1 3]
%     T = taskAgentMap==s;
%     
%     %if task label is 0 and score is 3, treat this as a negative
%     if s==3
%         T(negTasks, :) = -T(negTasks, :);
%     end
%     
%     %if task label is 1 and score is -1, treat this as a negative
%     if s==-1
%         T(posTasks, :) = -T(posTasks, :);
%     end
%     
%     %if score is 1, treat this as no co-occurrence
%     %     do nothing for score==1, skip this loop
%     
%     Adj = Adj + (T * T');
% end

%AN ALTERNATIVE MAY BE TO LOOK AT WHICH COMMUNITY THE AGENTS BELONG TO
%How to analyse if these communities are useful?
%Need to adjust the edge weights?

[I, J] = find(taskAgentMap);

% Tpos = sparse(I, J, reshape(Pi(1, 1, J)*Kappa(1)+Pi(2, 3, J)*Kappa(2), length(J), 1));
% Adj = Tpos*Tpos';
% Sim = abs(Adj);
% 
% totalsBySample = repmat(sum(Sim, 2), 1, nTasks);
% Sim = 2.* Sim ./ (totalsBySample + totalsBySample'); %multiply by 2 means maximum is 1
% 
% [Ps7, gS7] = commDetNMF(Sim, 100);
% mod = get_modularity(gS7, Sim, 'groups');
% display(num2str(mod));

Tpos = sparse(I, J, reshape(Pi(1, 1, J), length(J), 1));
Tpos = Tpos + sparse(I, J, reshape(Pi(2, 3, J), length(J), 1));
Tneg = sparse(I, J, reshape(Pi(1, 3, J), length(J), 1));
Tneg = Tneg + sparse(I, J, reshape(Pi(2, 1, J), length(J), 1));

% Adj = (Tpos-Tneg)*(Tpos-Tneg)';
Tpos = Tpos .^0.5;
Tneg = Tneg .^0.5;

Adj = Tpos*Tpos' - Tneg*Tneg';
Sim = abs(Adj);

% Sim = abs(Adj); % negative and positive connections both count, but the occurrence of both reduces the strength of the edges

totalsBySample = repmat(sum(Sim, 2), 1, nTasks);
Sim = 2.* Sim ./ (totalsBySample + totalsBySample'); %multiply by 2 means maximum is 1

[Ps5, gS5] = commDetNMF(Sim, 100);
mod = get_modularity(gS5, Sim, 'groups');
display(num2str(mod));
% 
% Adj = Tpos*Tpos';% + Tneg*Tneg';
% Sim = abs(Adj);
% 
% % Sim = abs(Adj); % negative and positive connections both count, but the occurrence of both reduces the strength of the edges
% 
% totalsBySample = repmat(sum(Sim, 2), 1, nTasks);
% Sim = 2.* Sim ./ (totalsBySample + totalsBySample'); %multiply by 2 means maximum is 1
% 
% [Ps6, gS6] = commDetNMF(Sim, 100);
% mod = get_modularity(gS6, Sim, 'groups');
% display(num2str(mod));
% 
% pause;

for c=1:numel(gS5)
    comm = gS5{c};
    display(Adj(comm, comm));
    display('press a key to continue');
    pause;
end

