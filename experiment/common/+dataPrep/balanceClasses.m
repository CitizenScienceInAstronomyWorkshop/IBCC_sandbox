function [idxsToKeep, subLabels] = balanceClasses(snBaseOutputs, labels, ...
    minAgentResp, minFreqAgents, keepUnlabelled, keepInfreqAgents, maxNoAssets, filterPos)
    %BALANCECLASSES 
    % Rebalance the class proportions. Assuming positive examples are far fewer
    % than negative ones, get rid of some negative examples.

    unknownAssets = find(labels==0);
    if nargin > 4 && keepUnlabelled==false
        unknownAssets = [];
    end

    if nargin < 6
        keepInfreqAgents = true;
    end

    posAssets = find(labels==2);
    posAssets = posAssets(ismember(posAssets, snBaseOutputs{2}));
    posAgents = sparse(snBaseOutputs{1}, 1, double(labels(snBaseOutputs{2})==2));

    if nargin > 7 && filterPos
        [idxsToKeep, subLabels, posAssets] = dataPrep.filterData(snBaseOutputs, posAssets, [],...
            posAgents, labels, minAgentResp, minFreqAgents, keepInfreqAgents, maxNoAssets);
    end
    
    display(['BalanceClasses: no. positive instances = ' num2str(length(posAssets))]);

    negAssets = find(labels==1);
    negAssets = negAssets(ismember(negAssets, snBaseOutputs{2}));
    negAssets = [negAssets; unknownAssets];
    if nargin > 6
        [idxsToKeep, subLabels, assetsToKeep] = dataPrep.filterData(snBaseOutputs, negAssets, posAssets,...
            posAgents, labels, minAgentResp, minFreqAgents, keepInfreqAgents, maxNoAssets);
    else
        [idxsToKeep, subLabels, assetsToKeep] = dataPrep.filterData(snBaseOutputs, negAssets, posAssets,...
        posAgents, labels, minAgentResp, minFreqAgents, keepInfreqAgents);
    end
    display(['BalanceClasses: no. negative instances = ' num2str(length(find(labels(assetsToKeep)==1)))]);
end

