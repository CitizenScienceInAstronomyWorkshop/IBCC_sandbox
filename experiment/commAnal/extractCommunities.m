function [g, P, avgPi, Alpha, Pi, observedAgents] = ...
    extractCommunities(Alpha, testData, snBaseOutputs, dynamic)

%Need indices for agents we observed.
observedAgents = unique(testData{1});%unique(testData{1});

observedTasks = unique(testData{2});%unique(testData{2});
nTasks = numel(observedTasks)

if (~exist('Alpha','var') || isempty(Alpha)) && ~isempty(agentRatings)
    Alpha = agentRatings{1}{3};    

%     missedSample = ismember(snBaseOutputs{1}, observedAgents) &...
%         ~ismember(snBaseOutputs{1}, testData{1}) & ...
%         labels(snBaseOutputs{2})~=0;    
%     
%     Alpha(labels(snBaseOutputs{2}(missedSample)), snBaseOutputs{3}(missedSample), snBaseOutputs{1}(missedSample)) = ...
%         Alpha(labels(snBaseOutputs{2}(missedSample)), snBaseOutputs{3}(missedSample), snBaseOutputs{1}(missedSample)) + 1;
end

%number of observations of each agent for positive examples and negative
%examples (supernova/not supernova)
% obsCountPos = sparse(1, testData{1}, double(labels(testData{2})==2) );
% obsCountNeg = sparse(1, testData{1}, double(labels(testData{2})==1) );


% %use Alpha already calculated.
if (nargin <4 || dynamic==false) && size(Alpha,3) > length(observedAgents)
    Alpha = Alpha(:, :, observedAgents); %check we haven't done this already!
end

nAgents = length(observedAgents);

if nargin > 3 && dynamic==true
    nAgents = size(Alpha, 3);
end

display(num2str(nAgents));

nClasses = size(Alpha,1);
nScores = size(Alpha,2);

normTerm = sum(Alpha, 2);
normTerm = repmat(normTerm, 1, nScores);
Pi = Alpha ./ normTerm;

%Need to know LnKappa
Kappa = zeros(1, nClasses) + 1./nClasses;
% negLabels = labels(observedTasks)==1;
% posLabels = labels(observedTasks)==2;
% Kappa(1) = 0.5;%sum(negLabels) / nTasks;
% Kappa(2) = 0.5;%sum(posLabels) / nTasks;

% minPi = reshape(repmat(min(Pi, 3), 1, nAgents), nClasses, nScores, nAgents);
% maxPi = reshape(repmat(max(Pi, 3), 1, nAgents), nClasses, nScores, nAgents);
% normPi = (Pi - minPi) ./ (maxPi - minPi);

% normTerm = repmat(sum(Pi, 2), [1, size(Pi,2), 1]);
% Pi = Pi ./ normTerm;

SquaredHD = squaredHellinger(nAgents, nClasses, nScores, Pi, Kappa);

% minSim = min(min(SquaredHD));
% maxSim = max(max(SquaredHD));
% normSim = (SquaredHD - minSim) ./ (maxSim - minSim);
% SquaredHD = normSim;

%try HD as dissimilarity? Or 1-similarity?
% Sim = exp(-SquaredHD*10);
Sim = exp(-SquaredHD*100);

% minSim = min(min(Sim));
% maxSim = max(max(Sim));
% normSim = (Sim - minSim) ./ (maxSim - minSim);
% Sim = normSim;
% 
% %co-occurrence normalisation
% totalsByAgent = sum(Sim, 2);
% totalsByAgent = repmat(totalsByAgent, 1, nAgents);
% Sim = Sim ./ (totalsByAgent + totalsByAgent');

maxRank = 3;
% W0 = zeros(size(Sim,1),maxRank);
% H0 = zeros(maxRank,size(Sim,2));
%P: soft membership score table. Rows->Agents. Columns->Communities.
%g: cell array containing lists of agents in each community
% [P, g] = commDetNMF(Sim, maxRank, W0, H0);
[P, g] = commDetNMF(Sim, maxRank);
display(g);
mod = get_modularity(g, Sim, 'groups');
display(['modularity: ' num2str(mod)]);
%number of positive examples observed by community
nObsPos = g;
nObsNeg = g;

scoreCounts = cell(1,numel(g));
scores = [-1 1 3; -1 1 3];

% for c=1:numel(g)
%     agentIdx = g{c};
%     
%     %chart the number of positive and negative examples in each community
%     nObsPos{c} = obsCountPos(observedAgents(agentIdx));
%     [N, X] = hist(nObsPos{c});
%     figure;
%     bar(X, N);
%     
%     nObsNeg{c} = obsCountNeg(observedAgents(agentIdx));
%     [N, X] = hist(nObsNeg{c});
%     figure;
%     bar(X, N);
    
%     %graph the counts of each score/true label combination for each
%     %community
%     scoreCounts{c} = [0 0 0; 0 0 0];
%     scoreCountAlt = [0 0 0; 0 0 0];
%     for j=1:size(scoreCounts{c},1)
%         for s=1:size(scoreCounts{c},2)
%             %indices in the set of classifications
%             dataIdxs = ismember(testData{1}, observedAgents(g{c}));
%             
%             %indices of classifications with score s
%             scoreMatches = testData{3}(dataIdxs)==scores(j,s);
%             
%             %the data sample idxs for selected classifications
%             commSamples = testData{2}(dataIdxs);
%             
%             %labels for the selected classifications
%             relevantLabels = labels(commSamples(scoreMatches));
%             
%             %count the relevant labels
%             scoreCounts{c}(j,s) = sum(relevantLabels==j);
%             
%             
%             for a=observedAgents(g{c})'
%                 dataIdxs = testData{1}==a;
%                 
%                 dataSamples = testData{2}(dataIdxs);
%                 relevantSamples = dataSamples(labels(dataSamples)==j);
%                 agentScores = testData{3}(dataIdxs);
%                 for n=relevantSamples'
%                     scoreIdxs = find(dataSamples==n);
%                     scoreCountAlt(j,s) = scoreCountAlt(j,s) + sum(agentScores(scoreIdxs)==scores(j,s));
%                 end
%             end
%             
%         end
%     end
%     
%     X = scores(1,:);
%     
%     scoreCountAlt
%     scoreCounts{c}
%     
% %     %normalise score counts
% %     scNormByClass = scoreCounts{c} ./ repmat(sum(scoreCounts{c}, 2), 1, 3);
% %     scNormByClass
% %     figure;
% %     bar(X, scNormByClass');
%     
%     scNormByScore = scoreCounts{c} ./ repmat(sum(scoreCounts{c}, 1), 2, 1);
%     scNormByScore
%     figure;
%     bar(X, scNormByScore');    
% end

avgPi = drawWeightedMeans(g, P, Pi);

end

% pause;
% 
% for c=1:numel(g)
%     commAgents = g{c};
%     
%     alphaComm = Alpha(:, :, commAgents);
%     
%     avgAlpha = sum(alphaComm, 3);
%     avgAlpha = avgAlpha ./ repmat(sum(avgAlpha, 2), 1, 3);
% 
%     figure;
%     X = [0 1];
%     bar3(X, avgAlpha);
% end
% 
% mod = get_modularity(g, Sim, 'groups');
% display(num2str(mod));



%-------------------------------------------------------------------------
%OTHER CLUSTERING METHODS

%number of negative examples observed by community

% Distance = 1-Sim;
% 
% Distance(Distance<0.00001) = 0;
% Distance(Distance>0.99999) = 1;

% Y = squareform(Distance);
% Z = linkage(Y, 'single');
% T = cluster(Z, 'maxclust', 5);
% 
% nComm = 1; nIt = 1;
% 
% distVec = squareform(Distance);
% 
% while nComm<5
%     
%     nIt
%     nIt = nIt + 1;
%     
%     [maxVal, maxIdx] = max(distVec);
%     distVec(maxIdx) = 1;
%     dist = squareform(distVec);
%     
%     comms = {[1]};
%     
%     for a=1:nAgents
%         
% %         a
%         
%         found = false;
%         for c=1:length(comms)
%             if ~isempty(find(comms{c}==a))
%                 found = true;
%                 break
%             end
%         end
%         if ~found
%             neighbours = find(dist(a,:)<1);
%             while ~isempty(neighbours) && ~found
%                 
%                 n = neighbours(1);
% %                 length(neighbours)
%                 for c=1:length(comms)
%                     if ~isempty(find(comms{c}==n))
%                         found = true;
%                         comms{c} = [comms{c} a];
%                         break
%                     end
%                 end
%                 
%                 neighbours = unique([neighbours(2:end) find(dist(n,:)<1)]);
%             end
%         end
%         if ~found
%             comms{nComm+1} = find(dist(a,:)<1);
%             nComm = length(comms);
%         end
%     end
%     nComm
% end
% 
% 
% 
