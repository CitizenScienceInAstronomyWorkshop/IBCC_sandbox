%Given a set of pi-communities, [P, g], compare the frequency of different
%samples in each pi-communities. Highlight those with high frequency and
%exclusivity to a particular pi-community

dataExpd = testData;

extraIdxs = ismember(snBaseOutputs{1}, testData{1}) & ~ismember(snBaseOutputs{2}, testData{2});

dataExpd{1} = [dataExpd{1}; snBaseOutputs{1}(extraIdxs)];
dataExpd{2} = [dataExpd{2}; snBaseOutputs{2}(extraIdxs)];
dataExpd{3} = [dataExpd{3}; snBaseOutputs{3}(extraIdxs)];

data = dataExpd;


[observedTasks sampleIdxToID sampleIDToIdx] = unique(data{2});
nTasks = numel(observedTasks);

[observedAgents, agentIdxToID, agentIDToIdx] = unique(data{1});
nAgents = length(observedAgents);

%For each pi-community, record the frequency of each sample
[p, testDataPiComm] = max(P(agentIDToIdx, :), [], 2);
gOverlay = sparse(sampleIDToIdx, testDataPiComm, ones(numel(data{2},1)));

gSizes = zeros(1, numel(g));
for c=1:numel(g)
    comm = g{c};
    gSizes(c) = numel(comm);
end

gOverlayCov = gOverlay ./ repmat(gSizes, size(gOverlay,1), 1);
gOverlayNorm = gOverlay ./ repmat(sum(gOverlay, 2), 1, numel(g));
gOverlayNorm = (gOverlayNorm ./ repmat(gSizes, size(gOverlay,1), 1)) .* sum(gSizes);

[maxVal, maxComm] = max(gOverlayNorm, [], 2);
[minVal, minComm] = min(gOverlayNorm, [], 2);
filteredIdxs = find(maxVal > 0.5 & sum(gOverlayCov,2)>0.06);% & maxComm ~= 3);

% figure;
% bar3(filteredIdxs, gOverlayNorm(filteredIdxs, :));

figure;
bar(1:numel(filteredIdxs), gOverlayNorm(filteredIdxs,:));

%coverage: how many of the community have seen this sample?
figure;
bar(1:numel(filteredIdxs), gOverlayCov(filteredIdxs,:));

figure;
[N,X] = hist(gOverlayNorm);
bar(X, N);


% gSizes = zeros(1, numel(g));
% for c=1:numel(g)
%     comm = g{c};
%     gSizes(c) = numel(comm);
% end
% 
% gOverlayNorm = gOverlay ./ repmat(gSizes, nTasks, 1);
% 
% [maxVal, maxComm] = max(gOverlayNorm, [], 2);
% [minVal, minComm] = min(gOverlayNorm, [], 2);
% filteredIdxs = find(maxVal > 0.027);
% 
% bar3(filteredIdxs, gOverlayNorm(filteredIdxs, :));
% figure;
% bar(1:numel(filteredIdxs), gOverlayNorm(filteredIdxs,:));