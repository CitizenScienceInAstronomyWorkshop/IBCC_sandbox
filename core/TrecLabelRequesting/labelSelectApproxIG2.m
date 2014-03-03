function samplesToRequest = labelSelectApproxIG2( nToTry, nRequests, X, P, T, labIdxs, nClasses )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    maxCert = max(P, [], 2);
    [leastCert, leastCertIdx] = sort(maxCert,1,'ascend');
    
    samplesToTry = leastCertIdx(1:nToTry);
    
    IG = zeros(size(P,1),nToTry);        
    for i=1:nToTry
        
        Pafter = zeros(size(P,1), size(P,2), nClasses);
        
        s = samplesToTry(i); %index to training samples
        P_s = P(s,:); %probability of the training sample class according to last run of classifier
               
 
        Xtrain = [X(labIdxs, :); X(s,:)];
        
        for j=1:nClasses
            newLabel = zeros(1, nClasses); newLabel(j) = 1;
            [P_class_t,P_feat_class_t] = multi_classifier(Xtrain, [T; newLabel] );
            pTgivC = multi_classify(P_class_t,P_feat_class_t,X);
            
            eAfter = -sum(pTgivC.*log(pTgivC),2);
            eAfter(isnan(eAfter)) = 0;
            
            eBefore = -sum(pTgivC.*log(P),2);
            eBefore(isnan(eBefore)) = 0;
            
            IGj = eBefore - eAfter;
            
            IG(:,i) = IG(:,i) + P_s(j).*IGj;
            Pafter(:,:,j) = pTgivC;
        end
        IG(s,i) = IG(s,i) - sum(sum(Pafter(s, :, :),3).*log(P_s));             
    end
    
    totalValue = zeros(1,nToTry);
    sampleSet = zeros(nRequests,nToTry);
    for i=1:length(samplesToTry)
        
        selected = i;
        
        for r=1:nRequests-1
            jointIGEst = sum(IG(samplesToTry,selected),2)./(r+1) + max(IG(samplesToTry,selected),[],2);
            [MInext MInextIdx] = min(jointIGEst);
            selected = [selected; MInextIdx];
        end
        
        jointIGEst = sum(IG(samplesToTry,selected),2)./(nRequests) + max(IG(samplesToTry,selected),[],2);
        sampleSet(:,i) = selected;
        totalValue(i) = sum(jointIGEst);
        
        display(['[' num2str(samplesToTry(sampleSet(:,i))') ']: ' num2str(sum(totalValue(i))) ', ' num2str(leastCert(i))]);        
    end
    
    [maxVal setToRequest] = max(totalValue);
    samplesToRequest = samplesToTry(sampleSet(:, setToRequest));
end

