function samplesToRequest = labelSelectApproxIG( nToTry, nRequests, P, currentTest, Alpha, currentLabels, nClasses )
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
               
        [AlphaMatList, lnAlphaMatList, AlphaMatTotalsList] = getAlphaInCmatForm(Alpha,currentTest+1);

        for j=1:nClasses
            Nu = sum(currentLabels,1)+1; Nu(j) = Nu(j) + 1;
            NuSum = sum(Nu);
            Kappa = Nu/NuSum;
            [IGj, eBefore, eAfter, Pafterj] = objectInfoGain(Kappa, P, P_s, currentTest+1, Alpha, ...
                j, currentTest(s,:)+1, AlphaMatList, lnAlphaMatList, AlphaMatTotalsList);
            IG(:,i) = IG(:,i) + P_s(j).*IGj;
            Pafter(:,:,j) = Pafterj;
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

