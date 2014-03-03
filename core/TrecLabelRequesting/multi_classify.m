function P_target = multi_classify(P_class,P_feat_class,X, useLikelihood, useNegFeatures)

if nargin < 4
    useLikelihood = false;
end

if nargin < 5
    useNegFeatures = true;
end

nclasses = length(P_class);

nSplits = 1;
splitSize = 500;
if size(X, 1) > splitSize
    nSplits = size(X, 1) ./ splitSize;
    nSplits = ceil(nSplits);
end

P_target = zeros(size(X,1), nclasses);

for s=0:nSplits-1
    
    if size(X,1)< s*splitSize+splitSize
        splitLimit = size(X,1) - s*splitSize;
    else
        splitLimit = splitSize;
    end

    Xsplit = X(s*splitSize+1:s*splitSize+splitLimit, :);
    
    logP = Xsplit*log(P_feat_class);
    if ~useLikelihood
        logP = logP +repmat(log(P_class),splitLimit,1);
    end
    if useNegFeatures
        logP = logP + (Xsplit==0)*log(1-P_feat_class);
    end

    logP = logP-repmat(max(logP,[],2),1,nclasses)+10;
 
    P_target_split = exp(logP)./repmat(sum(exp(logP),2),1,nclasses);
    
    P_target(s*splitSize+1: s*splitSize+splitLimit, :) = P_target_split;
end

