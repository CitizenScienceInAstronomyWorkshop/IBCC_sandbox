function [idxsToKeep, subLabels, assetsToKeep] = filterData(C, assetsToFilter, otherAssets,...
    agentValidity, labels, minAgentResp, minFreqAgents, keepInfreqAgents, maxNoAssets)
%FILTERDATA  
% Remove data points that only have responses from agents we don't have
% much data for.

if minFreqAgents==0 && minAgentResp==0
    assetsToKeep = assetsToFilter; %no filtering needed
else
    
    Vmat = sparse(C{1}, C{2}, agentValidity(C{1})>=minAgentResp );
    filterValidity = sum(Vmat(:,assetsToFilter),1);
    assetsToKeep = assetsToFilter( filterValidity>minFreqAgents );
%     assetsToKeep = [];
%     for n=assetsToFilter'
%         display(['filtering: ' num2str(n)]);
%         agents = unique(C{1}(C{2}==n));
%         agentCounts = agentValidity(agents);
%         if sum(agentCounts>=minAgentResp) > minFreqAgents
%             assetsToKeep = [assetsToKeep; n];
%         end
%     end
end

if nargin > 8 && maxNoAssets>0 && length(assetsToKeep)+length(otherAssets) > maxNoAssets
    start = 1;
    if maxNoAssets > length(assetsToKeep)
        maxNoAssets = length(assetsToKeep);
    end
    assetsToKeep = assetsToKeep(start:start+maxNoAssets-length(otherAssets));
end
assetsToKeep = [assetsToKeep; otherAssets];

%subsample labels; even if no filtering of assets takes place, only the labels
%for which at least one base classifier has provided a score will be included
subLabels = sparse(length(labels), 1);
subLabels(assetsToKeep) = labels(assetsToKeep);

%subsample decision data from base classifiers
idxsToKeep = ismember(C{2}, assetsToKeep);

if ~keepInfreqAgents
    idxsToKeep = idxsToKeep & (agentValidity(C{1})>minAgentResp);
end
end

