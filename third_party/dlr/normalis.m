%[Y] = normalise(X,D,mD,covD)
% Y : normalised data set
% X : input data set
% D : data set providing normal stats

function [Y] = normalise(X,D,mD,varD)

if (nargin < 4),
    c = var(D);
    m = mean(D);
else
    c = varD;
    m = mD;
end;

f = find(c==0);
c(f)=1;

X = X - ones(size(X,1),1)*m;	% removes mean
Y = X./(ones(size(X,1),1)*sqrt(c));	% unit variance

return;
