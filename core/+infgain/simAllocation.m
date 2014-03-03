%Goal: allocate 

A = 3; %number of allocations we'll actually calculate per user
H = 10; %horizon; number of allocations into the future we'll look; don't worry about probability of completion in this version

if ~exist('Alpha','var')
    nClasses = 2;
    N = 1000; %number of tasks
    K = 1000; %number of workers
    Alpha = zeros(nClasses, 3, K); %confusion matrices
    Nu = [10 10];
    ET = zeros(1, N); %probabilities of each task
else
    nClasses = size(Alpha,1);
    N = size(ET,2);
    K = size(Alpha,3);
end
Ass = zeros(N, K); %completed assignments of users to tasks

%%%%% BRUTE FORCE OPTIMAL METHOD UP TO H LOOK-AHEAD STEPS
%allocates every user available H times.

AllocUsers = zeros(1, h*K);
AllocTasks = zeros(1, h*K);

state = {Alpha, Nu, ET, Ass};

for h=1:A*K 
    IG = infoGain(1, state, H, 0, N, K, A); %the information gain of assigning task i to classifier k at lookahead step h

    [maxIGByUser, chosenTasks] = max(IG, [], 1);
    [maxIG, chosenUser] = max(maxIGByUser, [], 2);
    chosenTask = chosenTasks(chosenUser);

    AllocTasks(h) = chosenTask;
    AllocUsers(h) = chosenUser;
    state = updateState(state, chosenTask, chosenUser);
end

%%%%% PURE GREEDY (NO LOOK-AHEAD)

AllocUsers = zeros(1, h*K);
AllocTasks = zeros(1, h*K);
state = {Alpha, Nu, ET, Ass};
    
for h=1:A*K 
    [IG, newState, chosenTask, chosenUser] = greedyInfoGain(state, N, K, A); %the information gain of assigning task i to classifier k at lookahead step h

    AllocTasks(h) = chosenTask;
    AllocUsers(h) = chosenUser;
end

%%%%% ROLLOUT: OPTIMAL--ONE STEP, GREEDY--H STEPS

AllocUsers = zeros(1, h*K);
AllocTasks = zeros(1, h*K);

state = {Alpha, Nu, ET, Ass};

for h=1:A*K 
    IG = infoGain(1, state, 1, H, N, K, A); %the information gain of assigning task i to classifier k at lookahead step h

    [maxIGByUser, chosenTasks] = max(IG, [], 1);
    [maxIG, chosenUser] = max(maxIGByUser, [], 2);
    chosenTask = chosenTasks(chosenUser);

    AllocTasks(h) = chosenTask;
    AllocUsers(h) = chosenUser; 
    state = updateState(state, chosenTask, chosenUser);
end

%%%%% AUCTIONING USERS WITH ROLLOUT-BASED BIDS
%Does this come up with different answers? What is the advantage?
%Tries to avoid having to go so many levels deep. The long tail of tasks
%that have modest IG may require an expert whereas tasks where we can gain
%a lot could just be given to multiple novices. 

AllocUsers = zeros(1, h*K);
AllocTasks = zeros(1, h*K);

state = {Alpha, Nu, ET, Ass};

for h=1:A*K 
    IG = infoGain(1, state, 1, H, N, K, A); %the information gain of assigning task i to classifier k at lookahead step h

    sortedIG = sort(IG, 2);
    bids = (IG - repmat(sortedIG(:,2), 1, K)) * IG ./ repmat(sum(IG,2),1,K); %advantage over the second best assignment for that task
    
    [maxIGByUser, chosenTasks] = max(bids, [], 1);
    [maxIG, chosenUser] = max(maxIGByUser, [], 2);
    chosenTask = chosenTasks(chosenUser);

    AllocTasks(h) = chosenTask;
    AllocUsers(h) = chosenUser;
    state = updateState(state, chosenTask, chosenUser);
end

%%%% ITERATIVE UPDATES TO THE ALLOCATION?
