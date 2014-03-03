function samplesToRequest = labelSelectUncertGreedy( nRequests, X,...
    P, Alpha, Tcols, P_feat_class, P_class )
%labelSelectUncertDiff Select most uncertain but also try to select a set that is maximally
%different

    reqLimit = 5; %only ever ask for a document a maximum of 5 times
    
    samplesToRequest = zeros(1,nRequests);
    
    nClasses = size(P,2);
    nFeat = size(Alpha,3);
    featIdxs = 1:nFeat;
    
    for s=1:nRequests    
    
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

        %take the current most uncertain
        [chosenVal chosenIdx] = min(maxCert);    
        samplesToRequest(s) = unconfirmedIdxs(chosenIdx);
        
        %update the others with expected outcome from current choice
        features = X(samplesToRequest(s),:);
        
        updateIdxs = round(sub2ind([size(Alpha,2), size(Alpha,3)], features+1, featIdxs));
        for j=1:nClasses
            Alpha(j, updateIdxs) = Alpha(j, updateIdxs) + P(samplesToRequest(s), j);            

            alpha = reshape(Alpha(j, 2, :), [nFeat 1]);
            alphaSum = reshape(Alpha(j,1,:) + Alpha(j,2,:), [nFeat 1]);    
            P_feat_class(:,j) = alpha./alphaSum;
        end
        P = multi_classify(P_class,P_feat_class,X);
    end
end

