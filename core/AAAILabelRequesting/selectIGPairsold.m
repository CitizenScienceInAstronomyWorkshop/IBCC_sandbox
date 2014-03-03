function [ excludedWorkers, taskIdxs, permittedWorkers ] = selectIGPairs( ...
    method, K, nToTry, nToAllocate, X, nAllWorkers, P, combiner, excludedWorkers, ...
    Xfeat, ETnow, logJoint, centroidsFile)
%LABELSELECTIGPAIRS Select worker/object pairs
%   k - worker to allocate %nRequests - number of pairs to write
%   nToTry - number of possible pairs to test (we can't test all
%   possibilities so set this heuristically)
%   X - input data
%   P - probabilities of classes
%   combiner - the combiner object used to produce P. We need this to
%   simulate updates to decide which will be the most informative pair

%Method can be set to the folling strings:
% HFAS - hiring and firing with active selection
% AS - active selection but no hiring and firing
% OS - online screening of workers using a single confidence level, random
% doc selection
% Random - assigns docs at random with no hiring and firing
% For testing purposes:
% ASEntropy - quick active selection based on max entropy docs
% HFQuick - quick hiring and firing by looking at only one pair per worker

%For first n tasks, use HFAS so that all workers get some screening! So
%that people cannot say we need to do up-front screening only
if length(X{1}) < 50
    if strcmp(method,'Random') || strcmp(method,'ASEntropy')
        method = 'OS'; %perform screening but not AS
    elseif strcmp(method, 'AS') 
        method = 'HFAS';
    end
end

%add in the unknown worker 
K = [K; nAllWorkers+1];
nK = numel(K);

nObjs = size(P,1);
objsToTry = zeros(nK, nObjs); %record 1 for each worker/doc pair to test

Alpha_all = combiner.Alpha;
if iscell(Alpha_all)
    Alpha_all = Alpha_all{1};
end
nClasses = size(Alpha_all,1);
nScores = size(Alpha_all,2);

Alpha = zeros(nClasses, nScores, nK);
Pi = zeros(nClasses, nScores, nK);
lnPi = zeros(nClasses, nScores, nK);    

for n=1:nK
    k = K(n);

    possibleObjs = ones(1, size(P,1)); 
    prevResponses = (X{1}==k);

    if k==0 || k>nAllWorkers
        %combiner.Alpha0(:,:,1);
        Alpha(:,:,n) = combiner.Alpha0(:,:,1);%sum(combiner.Alpha{1},3) ./ size(combiner.Alpha{1},3);
    else
        lastIdx = find(prevResponses, 1, 'last');        
        Alpha(:,:,n) = Alpha_all(:,:,lastIdx);
    end
    
    Pi(:,:,n) = Alpha(:,:,n) ./ repmat(sum(Alpha(:,:,n),2), 1,nScores);
    lnPi(:,:,n) = psi(Alpha(:,:,n)) - psi(repmat(sum(Alpha(:,:,n),2), 1,nScores));

    possibleObjs(X{2}(prevResponses)) = 0;
    possibleObjs = find(possibleObjs);

    if nToTry>0
        poolSize = nToTry; %when we had a previous random sampling step, set pool 
        %size larger than nToTry and choose nToTry docs at random
        if poolSize>length(possibleObjs)
            poolSize = length(possibleObjs);
        end
    else
        poolSize = length(possibleObjs);
        nToTry = poolSize;
    end

    objsToTry_k = possibleObjs(1:poolSize);

    if length(objsToTry_k)>nToTry
        objsToTry_k = objsToTry_k(randi(poolSize, nToTry, 1));
    end

    objsToTry(n, objsToTry_k) = 1;
end

if strcmp(method,'AS') || strcmp(method,'HFAS') || strcmp(method, 'HFQuick')

    %perform clustering
    nClust = size(P,1)/25; %maximum number of clusters we want to deal with
    if nClust < nToAllocate
        nClust = nToAllocate;
    end

    % D = sparse(X{2}, X{1}, X{3}, size(Xfeat,2), nAllWorkers); %don't put this
    % in as it is reflected in P and has missing entries (gaps count as zeros
    % not as missing entries
    D = Xfeat';
    D = [D P];

    if exist(centroidsFile,'file')
        load(centroidsFile);
        if length(centroids)==nClust
            [clIdx, centroids, ~, ~] = kmeans(D, nClust, 'Start', centroids, 'EmptyAction', 'singleton');
        else
            [clIdx, centroids, ~, ~] = kmeans(D, nClust, 'EmptyAction', 'singleton');
        end
    else
        [clIdx, centroids, ~, ~] = kmeans(D, nClust, 'EmptyAction', 'singleton');
    end
    centroids = full(centroids);
    save(centroidsFile, 'centroids');

    clustIds = unique(clIdx);
    nClust = length(clustIds);
    chosenDocs = zeros(nK,nClust);
    %do sub-clustering of each cluster so that each worker is tested against
    %different subsets -> better coverage of input space esp. when choosing
    %between workers
    for cl=1:nClust
        members = find(clIdx==clustIds(cl));
    %     if length(members) >= nK
    %         [~, ~, ~, clDist] = kmeans(D(members, :), nK, 'EmptyAction', 'singleton');
    %         [~, chosenDocs(:, cl)] = min(clDist,[],1);
        memberIdx = randi(length(members),1,1);
        chosenDocs(:, cl) = members(memberIdx);
        memberIdx = 1;    
        while sum(objsToTry(:, chosenDocs(1, cl))==0)>0 && memberIdx<=length(members)
            chosenDocs(:, cl) = members(memberIdx);
            memberIdx = memberIdx+1;
        end

    %     else
    %         chosenDocs(:,cl) = members( randi(length(members), size(chosenDocs,1), 1) );
    %     end
    end
end

%Reasons for particular selections from the clusters:

%Use different subset for each agent: lots of similar agents -> no point
%repeating the calculation. By using different ones from each cluster, more
%chance of getting the good/bad ones from each cluster. So more average,
%less risky. Doesn't really increase diversity because we still select same
%number of decisions. Can increase optimality if we have choice between
%workers, i.e. not assigning one task per worker.
%equal chance of getting stuck on worse choices as finding better choices,
%so have to increase the number of choices by selecting between workers. 

%Reason to use same subset for each agent: if it is optimal subset, but the
%choice of subset is not very clear.

%Once a subset is chosen, it could be eliminated because: further
%assignments likely to cause duplication. This is not the case where one
%cluster has extremely high uncertainty or if we are doing silver tasking.

%The duplication risk is reduced if we select different members of each
%cluster and only a small number of concurrent assignments.

%Could reduce duplication risk by eliminating clusters in certain
%circumstances:
% When EIG over that cluster only is high (not silver tasking)
% Or when remaining average entropy per document remains high

%Alternatively, don't do any elimination because: there is always risk of
%sub-optimal choices for subsequent assignments even if not in same
%cluster; the most difficult areas need extra attention; small amounts of
%duplication is not a big waste; there's a chance that tasks are not
%completed, so make sure the important ones are done. 

%Can use hierarchical clustering to avoid duplication: choose a sub-cluster
%at random from within a larger cluster so overlap is reduced.

%Optimisation: if we are assigning each worker at least one task, and we 
%are eliminating clusters that are already assigned, do this first before
%extra tasks where we chose between workers.

IG = zeros(nK, length(objsToTry));

if strcmp(method,'AS') || strcmp(method,'HFAS') %complete IG table using all proposed documents from the clustering step
     for n=1:(nK-1)
        display(['evaluating agent ' num2str(K(n))]);
        k = K(n);
        
        IGn = zeros(1, length(chosenDocs(n,:)));
        Pin = Pi(:,:,n);
        parfor i=1:length(chosenDocs(n,:))
            if objsToTry(n, i)==0
                continue
            else
                display(num2str(i));
                IGn(i) = evaluatePair1(k, chosenDocs(n,i), P, logJoint, Pin, X, ETnow, combiner);
            end
        end
        
        IG(n,chosenDocs(n,:)) = IGn;
     end   

     priorHt = -sum(P(chosenDocs(1,:),:) .* log(P(chosenDocs(1,:),:)),2);
     [~, mostUncertain] = max(priorHt);    
     MU = chosenDocs(1,mostUncertain);
     display(['Unknown worker for most uncertain: ' num2str(MU)]);
     IG(nK,MU) = evaluatePair1(K(end), MU, P, logJoint, Pi(:,:,nK), X, ETnow, combiner);
    
elseif strcmp(method, 'HFQuick') %complete IG table with only one doc as a tester

    P(chosenDocs(1,1),:)

    priorHt = -sum(P(chosenDocs(1,:),:) .* log(P(chosenDocs(1,:),:)),2);
    [~, mostUncertain] = max(priorHt);
    MU = chosenDocs(1,mostUncertain);

    for n=1:(nK-1)
        display(['evaluating agent ' num2str(K(n))]);
        k = K(n);
        i = MU;
        if objsToTry(n, MU)==0
            i = find(objsToTry(n,:), 'first');
        end
%         display(num2str(i));
        IG(n, :) = evaluatePair1(k, i, P, logJoint, Pi(:,:,n), X, ETnow, combiner);
    end
    
elseif strcmp(method, 'ASEntropy')
    %complete IG table using maximum entropy
    IG = repmat(-sum(P .* log(P),2)', size(IG,1), 1);   
elseif strcmp(method, 'Random')
    IG = rand(size(IG));
    IG(nK,:) = 0; %don't allow the unkown worker to be assigned! No hiring and firing should take place.
elseif strcmp(method, 'OS')
    %Do online screening: assign each worker a confidence level, fill in IG matrix with this.
    %Don't differentiate between docs though.
    display('trust:');
    for n=1:nK
        trust = singleTrustValue(Alpha(:,:,n))
        IG(n,:) = trust;
    end
end

%select the optimal pairs of workers/documents.
%exclude those who are more useless than a random new member.
taskIdxs = zeros(1, nToAllocate);
permittedWorkers = zeros(1, nToAllocate);

remainIG = round(IG*100);

%ensure each worker has a task if better than the unknown worker
nExtraTasks = nToAllocate - (nK-1);
display(['n extra tasks ' num2str(nExtraTasks)]);

for n=1:nToAllocate
    
    [kMax kMaxIdx] = max(remainIG, [], 2);
    display(num2str(kMax'));
    [overallMax overallMaxIdx] = max(kMax);
    
    if overallMax==0 %unlikely case that we run out of tasks to allocate
        display('this allocation would give zero EIG');
        newIdx = overallMaxIdx;
        while ismember(newIdx, taskIdxs) 
            newIdx = newIdx +1;
            if  newIdx>size(IG,2)
                newIdx = 1;
            elseif newIdx==overallMaxIdx
                display('we allocated absolutely everything, so just going to allocate documents multiple times at random to new workers');
            end
        end
        overallMaxIdx = newIdx;
    end
    
    %compare with unknown worker    
    if strcmp(method,'OS')
        %need to choose document at random, we only do comparison to see if worker should be hired.
        %IG already set to the trust level (independent of document) above.
        availableI = find(IG(overallMaxIdx,:)>0);
        i = randi(numel(availableI), 1, 1);
        i = availableI(i);
    else
        i = kMaxIdx(overallMaxIdx);
    end
    
    if strcmp(method,'HFAS') || strcmp( method, 'HFQuick')
        %Methods that allow the unknown worker to be assigned instead of
        %the proposed worker (firing). OS method would already do this when
        %finding the max, since unknown worker is calculated in advance
        if IG(nK,i)==0%unknown worker has not been evaluated
            if strcmp( method, 'HFQuick') 
                i = MU;
            end
            IG(nK,i) = evaluatePair1(K(end), i, P, logJoint, Pi(:,:,nK), X, ETnow, combiner);
            remainIG(nK,i) = round(IG(nK,i)*100);

            Alpha(:,:,nK)
        end

        display(['unknown worker IG: ' num2str(remainIG(nK,i)) ' vs. ' num2str(remainIG(overallMaxIdx,i))]);

        if remainIG(nK,i) > overallMax
            overallMax = remainIG(nK,i);
            overallMaxIdx = nK;
            %chosen the unknown worker - don't give tasks to the others
            k = 0;
        elseif overallMaxIdx>nAllWorkers
            k = 0;
        else
            k = K(overallMaxIdx);
        end        
    else %method = AS, ASEntropy, Random; i.e. methods that don't fire workers online
        k = K(overallMaxIdx);
    end
        
    taskIdxs(n) = i;
    
    display([num2str(k) ': ' num2str(overallMax) ]);
    display(['EIG: ' num2str(overallMax) ', document ' num2str(i) ]);
         
    if ismember(k,permittedWorkers) || k==0
        nExtraTasks = nExtraTasks -1;
        display(['n extra tasks: ' num2str(nExtraTasks)]);
    end
        
    if k<=nAllWorkers && k > 0    
        %A way to favour workers that have nothing assigned to them and to
        %allocate to new workers when existing ones have enough tasks
        if nExtraTasks <= 0
            remainIG(overallMaxIdx,:) = 0;
        end
    end
    
    permittedWorkers(n) = k;
    
    remainIG(:,i) = 0; %remove the current object and look for new tasks
end

if size(K,1)>size(K,2)
    K = K';
end
for k=K
    if ~ismember(k, permittedWorkers) && k<=nAllWorkers
        display(['Excluding ' num2str(k)]);
        excludedWorkers = unique([excludedWorkers; k]);
    end
    
%     n = find(K==k);
%     Alpha(:,:,n)
end
 
if size(excludedWorkers,1) > size(excludedWorkers,2)
    excludedWorkers = excludedWorkers';
end

if size(taskIdxs,1) > size(taskIdxs,2)
    taskIdxs = taskIdxs';
end

if size(permittedWorkers,1) > size(permittedWorkers,2)
    permittedWorkers = permittedWorkers';
end

end

function trust=singleTrustValue(Alpha)
    nClasses = size(Alpha,1);

    nCorrect = sum(Alpha(sub2ind(size(Alpha), 1:nClasses, 1:nClasses)));
    total = sum(sum(Alpha));
    trust = nCorrect ./ total;
end

