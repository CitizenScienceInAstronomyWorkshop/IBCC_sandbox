function e = Tentropy(Alpha,Kappa,Cmat,PiIndx)
    %C can be scoreset or matrix
    nAssets = size(Cmat,1);
    nClasses = size(Alpha,1);
    nScores = size(Alpha,2);
    
    if size(Kappa,1)<nClasses
        Kappa = Kappa';
    end
    if size(Kappa,2)<nAssets
        Kappa = repmat(Kappa,1,nAssets);
    end
    
    nAgents = size(Alpha, 3);
    k = 1:nAgents;

    Pi = Alpha(:,:,k) ./ repmat(sum(Alpha(:,:,k),2), [1,size(Alpha,2),1]);

%     indx = sub2ind([nScores nAgents], Cvec, baseCl);     
    pCT = sparse([], [], [], nAssets, nClasses);
    
    for j=1:nClasses
        Pivals = Pi(j,PiIndx);
        PCgivT = reshape(Pivals,size(Cmat));
        pCT(:, j) = prod(PCgivT,2).* Kappa(j);        
% 
%         C4 = lnPi(j, PiIndx);
%         lnpCT(:,j) = sparse(objects, 1, C4, nAssets, 1) + log(Kappa(j));
    end
    pC = repmat(sum(pCT,2), 1, size(pCT,2));
    pTgivC = pCT./pC;
    
    e = -sum(sum(pTgivC(pTgivC~=0).*log(pTgivC(pTgivC~=0)),1),2);