function [IG, eBefore, eAfter, pTgivC] = objectInfoGain(Kappa,oldPT,newCTest,Alpha0,Ttrain,CTrain,AlphaMatList,lnAlphaMatList,AlphaMatTotalsList)
%The information gain for any objects on receiving updated true labels and base
%classifier outputs
%oldCTest - the data we want to evaluate entropy of labels given
%newCTest - the data we want to evaluate entropy of labels given, with
%updates
%CTrain - the data we want to calculate updated alpha parameters from on
%receiving a new label

%1. To find the IG caused by a new label or updated estimate of T, the parameters oldC should
% contain all base classifier outputs, while newC is empty. Ttrain 
% contains the new targets only, should be sparse and contain only new targets. LnKappa is the prior.
%2. To find the IG from a new base classifier output, oldC==[] and newC
%contains the new outputs, lnKappa is the probabilities given all previous
%data. newTargets is empty. oldET and Ttrain can be ignored.
%3. To update both, provide the new base classifier outputs in newC (only the new ones all previous ones), the
%old ones in oldC etc...

nClasses = size(Alpha0,1);
nScores = size(Alpha0,2);
nAgents = size(Alpha0,3);
nAssets = size(newCTest,1);

[lnAlphaMatList, AlphaMatTotalsList] = updateAlphaInCmatForm(AlphaMatList, lnAlphaMatList, AlphaMatTotalsList, CTrain, Ttrain,newCTest);



PiMatList = cell(1,nClasses);

pCT = zeros(nAssets, nClasses);
for j=1:nClasses
    if length(AlphaMatTotalsList{j})>1
        AlphaTotals = sum(log(AlphaMatTotalsList{j}),2);%prod(AlphaMatTotalsList{j}, 2);
    else
        AlphaTotals = nAgents .* log(AlphaMatTotalsList{j});
    end
    PiMatList{j} = sum(lnAlphaMatList{j},2);
    pCT(:,j) = exp(PiMatList{j} - AlphaTotals + log(Kappa(j)) );
end
pC = repmat(sum(pCT,2), 1, size(pCT,2));
pTgivC = pCT./pC;

eAfter = -sum(pTgivC.*log(pTgivC),2);
eAfter(isnan(eAfter)) = 0;
% AlphaUpdated = Alpha0; 
% AlphaUpdated(Ttrain, sub2ind([size(Alpha0,2), size(Alpha0,3)], CTrain, 1:nAgents)) ...
%     = AlphaUpdated(Ttrain, sub2ind([size(Alpha0,2), size(Alpha0,3)], CTrain, 1:nAgents)) + 1;
% eAlphaAfter = sum(sum(dirEntropy(AlphaUpdated),1),3);
% eAfter = sum(eAfter);
% 
% if nargin>5 && ~isempty(Ttrain)
%     Count = zeros(nClasses, nScores, nAgents);     
%     if ~iscell(CTrain)
%         for j=1:nClasses
%             for s=1:nScores
%                 Count(j, s, :) = sum(repmat(Ttrain(:,j)', nAgents, 1) .* CTrain'==s, 2);
%             end
%         end
%     else
%         for j=1:nClasses
%             for s=1:nScores
%                 Count(j, s, :) = sparse(CTrain{1}, ones(length(CTrain{1})), ...
%                     repmat(Ttrain(j,CTrain{2}), obj.nAgents, 1) .* CTrain{3}==s, ...
%                     nClasses, nScores);
%             end
%         end        
%     end
%     
%     newAlpha = Alpha0 + Count;
% end

% Cvec = oldCTest{3};
% baseCl = oldCTest{1};
% objects = oldCTest{2};
% eBefore = Tentropy(oldAlpha,lnKappa,objects,baseCl,Cvec);

% [lnAlphaMatList, AlphaMatTotalsList] = updateAlphaInCmatForm(AlphaMatList, lnAlphaMatList, AlphaMatTotalsList, CTrain, oldPT,newCTest);
% PiMatList = cell(1,nClasses);
% 
% pCT = zeros(nAssets, nClasses);
% for j=1:nClasses
%     if length(AlphaMatTotalsList{j})>1
%         AlphaTotals = sum(log(AlphaMatTotalsList{j}),2);%prod(AlphaMatTotalsList{j}, 2);
%     else
%         AlphaTotals = nAgents .* log(AlphaMatTotalsList{j});
%     end
%     PiMatList{j} = sum(lnAlphaMatList{j},2);
%     pCT(:,j) = exp(PiMatList{j} - AlphaTotals) .* Kappa(j);
% end
% pC = repmat(sum(pCT,2), 1, size(pCT,2));
% pTgivC = pCT./pC;
% eBefore = -sum(sum(pTgivC(pTgivC~=0).*log(pTgivC(pTgivC~=0)),1),2);

% [lnAlphaMatList, AlphaMatTotalsList] = updateAlphaInCmatForm(AlphaMatList, lnAlphaMatList, AlphaMatTotalsList, CTrain, PTtrain,newCTest);
% PiMatList = cell(1,nClasses);
% 
% pCTPrior = zeros(nAssets, nClasses);
% for j=1:nClasses
%     if length(AlphaMatTotalsList{j})>1
%         AlphaTotals = sum(log(AlphaMatTotalsList{j}),2);%prod(AlphaMatTotalsList{j}, 2);
%     else
%         AlphaTotals = nAgents .* log(AlphaMatTotalsList{j});
%     end
%     PiMatList{j} = sum(lnAlphaMatList{j},2);
%     pCTPrior(:,j) = exp(PiMatList{j} - AlphaTotals + log(Kappa(j)) );
% end
% pCPrior = repmat(sum(pCTPrior,2), 1, size(pCTPrior,2));
% pTgivCPrior = pCTPrior./pCPrior;
eBefore = -sum(pTgivC.*log(oldPT),2);
eBefore(isnan(eBefore)) = 0;

% %this is probably not quite right...
% AlphaBefore = Alpha0; 
% AlphaBefore(:, sub2ind([size(Alpha0,2), size(Alpha0,3)], CTrain, 1:nAgents)) ...
%      = AlphaBefore(:, sub2ind([size(Alpha0,2), size(Alpha0,3)], CTrain, 1:nAgents)) + repmat(PTtrain', 1, nAgents);
%  
% eAlphaBefore = sum(sum(dirEntropy(AlphaBefore),1),3);
% eBefore = sum(eBefore);


% Cvec = newCTest{3};
% baseCl = newCTest{1};
% objects = newCTest{2};
% eAfter = Tentropy(newAlpha,Kappa,newCTest,PiIndx);
% eAlphaBefore - eAlphaAfter
IG = eBefore - eAfter;


end
