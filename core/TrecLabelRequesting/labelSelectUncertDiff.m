function samplesToRequest = labelSelectUncertDiff( nToTry, nRequests, X, P, Alpha, Tcols )
%labelSelectUncertDiff Select most uncertain but also try to select a set that is maximally
%different

    reqLimit = 3; %only ever ask for a document a maximum of 5 times

    pLabCorrect = 0.9; pIncorrect = 1-pLabCorrect;
    T = sparse( Tcols(:,2), Tcols(:,3), 1, size(P,1), size(P,2) );
    labelCounts = repmat(sum(T,2), 1, size(T,2));
    T(T==0) = labelCounts(T==0) .* log(pIncorrect);
    T(T>0) = labelCounts(T>0) .* log(pLabCorrect);
    P = exp(log(P) + T); 
    P = P ./ repmat( sum(P,2), 1, size(P,2) ); %incorporate labels
    maxCert = max(P, [], 2);
        
    
    %sometimes uncertainty collapses to zero despite seeing only a few
    %labels. At this point we just want to pick at random from different
    %clusters of documents. However, avoid picking the same documents again
    %as they have led us into this sticky mess and seeing them again won't
    %help at this stage.
    dontAskIdxs = [Tcols( maxCert(Tcols(:,2))==1, 2); find(labelCounts(:,1)>reqLimit)];
    unconfirmedIdxs = setdiff( 1:size(P,1), dontAskIdxs )';
    maxCert = maxCert(unconfirmedIdxs);
    
    
    %There seems to be a trade-off between exploring new parts of the
    %feature space and correcting previous errors from earlier dodgy labels
    %before they compound themselves. Adding this part in tries to explore
    %features that don't yet have many known labels but seems to lead to
    %poorer results.
%     countTotals = sum(repmat(labelCounts(:,1), 1, size(X,2)).*X, 1);%sum(X(labelCounts(:,1)>0,:),1);
%     possibleCounts = sum(X,1);
%     countFraction = countTotals ./ possibleCounts;
%     docCountTotals = X*countFraction'; %documents with features that already have high counts should be weighted less
%     maxCert = maxCert ./ docCountTotals(unconfirmedIdxs);
    
    [leastCert, leastCertIdx] = sort(maxCert,1,'ascend');    
    if nToTry>length(maxCert)
        nToTry = length(maxCert);
    end
    samplesToTry = leastCertIdx(1:nToTry);
    
%     chosenIdxs = 1:nRequests;
    
    [clusterIdx, centroids, sumDist, dist] = kmeans( ...
        X(samplesToTry,:), nRequests, 'EmptyAction', 'singleton');

    chosenIdxs = zeros(1, nRequests);
    for k=1:nRequests
        idxs = find(clusterIdx==k);
        [tmpVal chosenIdx] = min(leastCert(idxs));
        chosenIdxs(k) = idxs(chosenIdx);
    end
    
    %convert back to global doc index
    samplesToRequest = unconfirmedIdxs(samplesToTry(chosenIdxs));
end

