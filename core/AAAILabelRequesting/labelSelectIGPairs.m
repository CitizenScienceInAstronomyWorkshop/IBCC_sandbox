function [ excludedWorkers, taskIdxs, permittedWorkers ] = labelSelectIGPairs( ...
    k,nToTry, X, nAllWorkers, P, nClasses, combiner, excludedWorkers, Tvec, bccSettings, nFeat, featSettings, ...
            Xfeat )
%LABELSELECTIGPAIRS Select worker/object pairs
%   k - worker to allocate %nRequests - number of pairs to write
%   nToTry - number of possible pairs to test (we can't test all
%   possibilities so set this heuristically)
%   X - input data
%   P - probabilities of classes
%   combiner - the combiner object used to produce P. We need this to
%   simulate updates to decide which will be the most informative pair
%   workerIdx - indexes into columns of X that correspond to workers rather
%   than features.

%Entopy calculation:
% p(t|c) and log(p(t|c)) are given by the variational method. 
% The problem: if all base classifiers are initially thought to be utterly
% useless then we never predict an information gain over the labels. We
% would predict an information gain over pi, however. TBH this situation is
% nonsensical as we always need to provide some known labels or informative
% priors to get a result, so this situation won't occur in practice.
logP = log(P);
logP(isinf(logP)) = 0;
Hnow_t = -sum(P .* logP, 2);
ETnow = combiner.combinedPost;

%single worker allocation
IG = zeros(size(P,1), 1);

%Batch
% workers = unique(X{1}(end-nRequests+1:end));
% nWorkers = length(workers);
% IG = zeros(size(P,1), nWorkers+1);

% for i_k=1:(nWorkers+1)
%     if i_k > nWorkers
%         k = nAllWorkers + 1;
%     else
%         k = workers(i_k);
%     end
 
if k==0
    k = nAllWorkers+1;
end

if ismember(k, excludedWorkers) 
    display(['Worker ' num2str(k) ' has been excluded.']);
    taskIdxs = [];
    permittedWorkers = [];
    return;
end
if nToTry==0
    taskIdxs = [];
    permittedWorkers = [];
    return;
end
possibleObjs = ones(1, size(P,1)); 
prevResponses = (X{1}==k);

if sum(prevResponses)<10 && k <= nAllWorkers
    taskIdxs = [];
    permittedWorkers = [];
    return
end
% %speed hack
% display('hack!')
Hnow_poss = Hnow_t(possibleObjs);
[ent, maxEntIdx] = sort(Hnow_poss, 'descend');
% endIdx = find(ent==ent(100),1,'last');
% taskIdxs = maxEntIdx(randi(endIdx, 100, 1));

if k>nAllWorkers %k>nWorkers
    Alpha = combiner.Alpha0(:,:,1);
else
    lastIdx = find(prevResponses, 1, 'last');        
    Alpha = combiner.Alpha;
    if iscell(Alpha)
        Alpha = Alpha{1};
    end
    Alpha = Alpha(:,:,lastIdx);
end
nScores = size(Alpha,2);
Pi = Alpha ./ repmat(sum(Alpha,2), 1,nScores);

% if length(prevResponses)>10 && max(max(Pi))<= 11/111
%     excludedWorkers = [excludedWorkers k];
%     display('exclusion hack');
%     permittedWorkers = [];
% else
%     permittedWorkers = k;
% end
% return;
%end of hack

possibleObjs(X{2}(prevResponses)) = 0;
possibleObjs = find(possibleObjs);
nPoss = length(possibleObjs);
Hnow_poss = Hnow_t(possibleObjs);
[~, maxEntIdx] = sort(Hnow_poss, 'descend');

poolSize = nToTry * 10;
if poolSize>length(maxEntIdx)
    poolSize = length(maxEntIdx);
end
objsToTry = maxEntIdx(1:poolSize);

if length(objsToTry)>nToTry
    objsToTry = objsToTry(randi(length(objsToTry), 1, nToTry));
end

% display(['*** testing worker ' num2str(k)]);

for i=possibleObjs(objsToTry)
    IG(i) = evaluatePair1(k,i,P, Pi, bccSettings, nFeat, featSettings, ...
                Xfeat, nAllWorkers, Tvec, X);
end
% end

%select the optimal object for each worker
[objIG objIdxs] = max(IG, [], 1);

%exclude those who are more useless than a random new member. Multiply by
%0.9 to make sure that it isn't a random effect of the choice of different
%samples to test

thetaI = zeros(size(objIG));
    
Pi_k = Pi;
Alpha_k = Alpha;

Alpha = combiner.Alpha0(:,:,1);

mk = max(max(Alpha_k));
mu = max(max(Alpha));

Pi = Alpha ./ repmat(sum(Alpha,2), 1,nScores);

%test the known agents against the unknown - are they worth keeping?
for i_i=1:length(objIdxs-1) %skip the last one, it's the unknown agent anyway
    i = objIdxs(i_i);

    Pti = P(i,:);

    %simulate the updates     
    Hnew = zeros(1,nClasses);
    HnewPc = zeros(1,nClasses);

    pC = sum(Pi .* repmat(Pti', 1, size(Pi,2)), 1);
    if pC(1)~=pC(2) %skip this if they're all the same. A fair comparison is with the same response as the worker being compared.
        [~, subsetC] = sort(pC, 'descend');
        subsetC = subsetC(1); 
    end
    for c=subsetC
        X{1} = [X{1}; nAllWorkers+1];
        X{2} = [X{2}; i];
        X{3} = [X{3}; c];
        combCopy = combiners.bcc.MixedIbccVb(bccSettings, nFeat, featSettings, ...
            Xfeat, nAllWorkers+1, 1, Tvec, [], nClasses, nClasses);
        combCopy.initAll = true;
        combCopy.combinedPost = ETnow;
        combCopy.combineDecisions(X); %check that we don't re-initialise ET       
        Pnew = exp(combCopy.logJoint);
        logPnew = combCopy.logJoint;
        logPnew(isinf(logPnew)) = 0;

        Hnew(c) = -sum(sum(Pnew .* logPnew, 2));
        Hnow = -sum(sum(Pnew .* logP', 2));        
        HnewPc(c) = (Hnow-Hnew(c)) .*pC(c);

        X{1} = X{1}(1:(end-1));
        X{2} = X{2}(1:(end-1));
        X{3} = X{3}(1:(end-1));
    end
    
    thetaI(i_i) = sum(HnewPc) ./ sum(pC(subsetC));
end

%batch method
% taskIdxs = unique(objIdxs);
% taskIdxs = taskIdxs';
% permittedWorkers = cell(length(taskIdxs), 1);
% workers(nWorkers+1) = nAllWorkers+1;
% for t=1:length(taskIdxs)
%     workers(objIdxs==taskIdxs(t))
%     permittedWorkers{t} = workers(objIdxs==taskIdxs(t));
% end

% igCols = dlmread('/homes/49/edwin/data/aaai/test_hiring/workerIG.csv');
% igCols = [igCols; [k objIG thetaI*0.999]];
% dlmwrite('/homes/49/edwin/data/aaai/test_hiring/workerIG.csv', igCols);

% if sum(prevResponses)<10
%     taskIdxs = [];
%     permittedWorkers = [];
%     return
% end

%possibly expect a higher standard from approved members? But how to determine it?
%use a threshold of 0.95 to allow for the approximate nature of the result
%from VB
if thetaI > 0.001
    thresholdedTasks = objIG < (thetaI*0.999);
else
    thresholdedTasks = 0;
end
display([num2str(k) ': ' num2str(objIG) ' vs. ' num2str(thetaI) ' ( < ' num2str(thetaI*0.999) ')']);
%ban these workers who have now dropped below usefulness
%Batch mode
% excludedWorkers = unique([excludedWorkers; workers(thresholdedTasks==0)]);

%single worker mode
if thresholdedTasks && sum(prevResponses)>=10
    display(['Excluding ' num2str(k)]);
    excludedWorkers = unique([excludedWorkers; k]);
    taskIdxs = [];
    permittedWorkers = [];
else
    %single worker method
    taskIdxs = objIdxs;
    permittedWorkers = k;
end

%TODO! avoid giving the same object to several workers at once
%Possible procedures:
% 1. a) give duplicate top task to whoever starts first. This avoids waiting
% for a long time for someone to do the task, avoids duplicating the task.
% However, may cause the best worker not to get the best task.
% b) give top task only to the highest IG worker. Give the others the next
% best task. Could block certain tasks if the top worker is not taking the
% task; not optimal because choice of who to assign second best task etc.
% is arbitrary.
% 2. a) Submit alternative tasks as soon as a task is taken OR
% b) submit at least two possible tasks for each worker if there is a clash.
% c) just wait till that task is complete. This could result in blocking
% some workers if they have to wait for a task to complete.
% On average, choosing option 1(a) and 2(a) or (b) is as good as choosing
% 1(b) then 2(a) or (b) due to greedy algorithm. However, 1(a) prevents
% task blocking.
%Choosing 2(c) leads closer to optimal allocation, but in practice it may
%discourage good workers who are starved of tasks so lead to suboptimal
%take-up of tasks. This depends on how long it takes before the loop completes
%and a new task is posted.
% 2(b) is computationally efficient and allows workers some choice of tasks (motivation?).
% 2(a) ensures the top task is taken but requires tasks to be updated at an
% additional step when top task is taken - so possibly not much advantage
% over a. 
%Go with 1(a) and 2(c) for now as it is simplest to get running - other 
% options require complex analysis. With 2(b) additional tasks in event of 
% clash are determined as follows: while number of tasks on offer to a 
% worker < minConcurrent and IG of proposed pair > lowest first choice 
% pair: if there is a clash, add next best task. 
%This just ensures that all workers have an equal chance to take on tasks
%in case running this loop takes a long time, allowing some concurrency for speed. 
%Can't have too much concurrency otherwise we move away from our improved
%expected information gain. If we don't care about the runtime of this
%script discouraging workers


end

