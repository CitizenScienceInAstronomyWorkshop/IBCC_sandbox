function [lnAlphaMatList, AlphaMatTotalsList] = updateAlphaInCmatForm(AlphaMatList,lnAlphaMatList,AlphaMatTotalsList,newCrow,T,Cmat)

newCmat = repmat(newCrow, size(Cmat,1), 1);

if length(T)==1
    j = T;
    lnAlphaMatList{j}(Cmat==newCmat) = log(AlphaMatList{j}(Cmat==newCmat)+1);
    AlphaMatTotalsList{j} = AlphaMatTotalsList{j}+1;
else
    for j=1:length(AlphaMatList)
        lnAlphaMatList{j}(Cmat==newCmat) = log(AlphaMatList{j}(Cmat==newCmat)+T(j));
        AlphaMatTotalsList{j} = AlphaMatTotalsList{j}+T(j);
    end
end