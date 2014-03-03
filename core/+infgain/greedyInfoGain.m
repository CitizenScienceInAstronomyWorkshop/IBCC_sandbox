function [chosenTask, chosenUser, maxIG, newState] = greedyInfoGain(state, N, K, A, nToSelect, updateState)
    if nargin < 4
        nToSelect = 1;
    end
        
    IG = zeros(N, K);
    
    Ass = state{4};
    ET = state{3};
    
    for i=1:N
        
        if ~isempty(find(ET(:,i)==1, 1))
            continue;
        end
        
        ig = batchInfoGain(i, state, true, true); %calculate here using state
        IG(i,:) = ig;
        IG(i, Ass(i,:)>=A) = 0;            
    end
    
    [maxIGByUser, chosenTasks] = max(IG, [], 1);
    [maxIGByUser, uId] = sort(maxIGByUser);
    chosenTasks = chosenTasks(uId);
    chosenTask = chosenTasks(1:nToSelect);
    chosenUser = uId(1:nToSelect);
    maxIG = maxIGByUser(1:nToSelect);
    
    if nargin > 5 && updateState
        newState = updateState(state, chosenTask, chosenUser);%state after making this allocation
    end
end
