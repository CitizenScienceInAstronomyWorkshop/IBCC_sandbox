function [AlphaMatList, lnAlphaMatList, AlphaMatTotalsList] = getAlphaInCmatForm(Alpha,Cmat)

nScores = size(Alpha,2);
nAgents = size(Alpha,3);
nClasses = size(Alpha,1);

AlphaMatList = cell(1,nClasses);
lnAlphaMatList = cell(1,nClasses);
AlphaMatTotalsList = cell(1,nClasses);

for j=1:nClasses
    AlphaMat = Cmat;
    totalSCount = zeros(nScores,nAgents);

    for s=1:nScores
        totalSCount(s,:) = reshape(Alpha(j,s,:), 1, nAgents);
        totalSCountRep = repmat(totalSCount(s,:), size(Cmat,1), 1);
        AlphaMat(Cmat==s-1) = totalSCountRep(Cmat==s-1);
    end

    totals = sum(totalSCount,1);
    if totals==totals(1)
        AlphaMatTotals = totals(1);
    else
        AlphaMatTotals = repmat(totals, size(Cmat,1), 1);
    end
    
    AlphaMatList{j} = AlphaMat;
    lnAlphaMatList{j} = log(AlphaMat);
    AlphaMatTotalsList{j} = AlphaMatTotals;
end
end