function samplesToRequest = labelSelectUncert( nRequests, X,...
    P, Alpha, Tcols )
%labelSelectUncertDiff Select most uncertain but also try to select a set that is maximally
%different

    reqLimit = 5; %only ever ask for a document a maximum of 5 times
    
    samplesToRequest = zeros(1,nRequests);   

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

    for s=1:nRequests
        %take the current most uncertain
        [chosenVal chosenIdx] = min(maxCert);    
        samplesToRequest(s) = unconfirmedIdxs(chosenIdx);
        maxCert(chosenIdx) = 1; %make sure we don't choose the same one twice
    end
end

