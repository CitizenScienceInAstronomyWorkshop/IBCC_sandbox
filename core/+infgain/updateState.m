function newState = updateState(oldState, i, k)
%UPDATESTATE Updates the distribution over the state variables after this
%task assignment has been completed, without knowing the results. 
%i is the chosen task; k is the chosen classifier.
    Ass = oldState{4};
    Ass(i,k) = Ass(i,k) + 1;
    
    
    Alpha = oldState{1};
    Nu = oldState{2};
    ET = oldState{3};
    
    pCT = pClassifierOutput(oldState, i,k);
    Alpha(:,:,k) = Alpha(:,:,k) + pCT;
    
    pCT = pClassifierOutput(oldState, i, k);
    pC = repmat(sum(pCT, 1), size(Alpha,1), 1);
    
    Nu = Nu - ET(:,i)';
    ET(:,i) = sum(pCT ./ pC, 2);
    Nu = Nu + ET(:,i)';
    
    newState = {Alpha, Nu, ET, Ass};    
end

