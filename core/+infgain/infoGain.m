function IG = infoGain(h, state, Hbranch, Hgreedy, N, K, A)
    IG = zeros(N, K);
        
    for i=1:N
        for k=1:K
            
            Ass = state{4};
            if Ass(i,k)>=A
                continue;
            end
            
            IG(i,k) = immediateInfoGain(i,k,state); %calculate here using state
            newState = updateState(state, i, k);%state after making this allocation
            if h < Hbranch
                IG_next = infoGain(h+1, newState, H, N, K, A);
                
                [maxIGByUser, chosenTasks] = max(IG_next, [], 1);
                [maxIG, chosenUser] = max(maxIGByUser, [], 2);
                
                %Assume we chose the maximum IG path
                IG(i,k) = IG(i,k) + maxIG;
            else
                for h=2:Hgreedy
                    [IG_next, newState, chosenTask, chosenUser] = greedyInfoGain(newState, N, K, A); %the information gain of assigning task i to classifier k at lookahead step h
                    IG(i,k) = IG(i,k) + IG_next;
                end
            end
        end
    end
end
