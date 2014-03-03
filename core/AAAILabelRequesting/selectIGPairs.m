function [ excludedWorkers, taskIdxs, permittedWorkers, Alpha, agentScores, unknownScores ] = selectIGPairs( ...
    method, K, nToTry, nToAllocate, X, nAllWorkers, P, combiner, excludedWorkers, ...
    Xfeat, ETnow, logJoint, centroidsFile, IGTargetsOnly)
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

if ~exist('IGTargetsOnly','var')
    IGTargetsOnly = true;
end

%For first n tasks, use HFAS so that all workers get some screening! So
%that people cannot say we need to do up-front screening only
%  if length(X{1}) < 50
%      if strcmp(method,'Random') ||strcmp(method,'Random2') || strcmp(method,'ASEntropy')
%          display('Random should choose any random agent, not from a hired pool');
%  %          method = 'OS'; %perform screening but not AS
%  %  	display('Should we do some pre-screening with real workers for the random selection method? Pre-screening using gold labels could give an unfair advantage unless gold are available to all. Can apply pre-screening to random method after collecting labels (i.e. on same collected data), so do not bother with it now');
%      elseif strcmp(method, 'AS') 
%          display('AS should choose any random agent, not from a hired pool');
%  %  	display('Should we do some pre-screening with real workers for the Active Selection method? Pre-screening using gold labels could give an unfair advantage unless gold are available to all. Can possibly apply pre-screening to AS method after collecting labels.');
%  %          method = 'HFAS';
%      else
%  %  	display('How do we compare with the pre-screening method esp. with real workers? Do the pre-screening after label collection for ALL methods, i.e. on the same set of labels. Not quite fair for HF and AS though. Do a second run with pre-screening using documents NOT from the main document set.');
%      end
%  end

%  display('Forget pre-screening now as this could equally be applied to all methods. Look to ensure that advantage of HFFullyDyn is not purely because it is an alternative to pre-screening; show the agents changing in reality');

%add in the unknown worker 
K(K==0) = nAllWorkers+1;
K = [K; nAllWorkers+1];
nK = numel(K);

nObjs = size(P,1);
objsToTry = ones(nK, nObjs); %record 1 for each worker/doc pair to test

Nu = combiner.Nu;
AlphaFeat = combiner.featAlpha;

Alpha_all = combiner.Alpha;
AlphaSum_all = combiner.AlphaSum;
if iscell(Alpha_all)
    Alpha_all = Alpha_all{1};
end
nClasses = size(Alpha_all,1);
nScores = size(Alpha_all,2);

Alpha = zeros(nClasses, nScores, nK);
AlphaSum = zeros(nClasses, nScores, nK);
Pi = zeros(nClasses, nScores, nK);
lnPi = zeros(nClasses, nScores, nK);    

for n=1:nK
    k = K(n);

%     possibleObjs = ones(1, size(P,1)); 
    prevResponses = (X{1}==k);

    if k==0 || k>nAllWorkers
        %combiner.Alpha0(:,:,1);
        Alpha(:,:,n) = combiner.Alpha0(:,:,1);%sum(combiner.Alpha{1},3) ./ size(combiner.Alpha{1},3);
        AlphaSum(:,:,n) = repmat(sum(Alpha(:,:,n),2), 1, size(Alpha,2));
    else
        lastIdx = find(prevResponses, 1, 'last');        
        Alpha(:,:,n) = Alpha_all(:,:,lastIdx);
        AlphaSum(:,:,n) = AlphaSum_all(:,:,lastIdx);
    end
    
    Pi(:,:,n) = Alpha(:,:,n) ./ AlphaSum(:,:,n);%repmat(sum(Alpha(:,:,n),2), 1,nScores);
    lnPi(:,:,n) = psi(Alpha(:,:,n)) - psi(AlphaSum(:,:,n));%psi(repmat(sum(Alpha(:,:,n),2), 1,nScores));

%     possibleObjs(X{2}(prevResponses)) = logical(false);
%     possibleObjs = find(possibleObjs);

%     if nToTry>0
%         poolSize = nToTry; %when we had a previous random sampling step, set pool 
        %size larger than nToTry and choose nToTry docs at random
%         if poolSize>length(possibleObjs)
%             poolSize = length(possibleObjs);
%         end
%     else
%         poolSize = size(P,1);%length(possibleObjs);
%         nToTry = poolSize;
%     end

%     objsToTry_k = possibleObjs(1:poolSize);

%     if length(objsToTry_k)>nToTry
%         objsToTry_k = objsToTry_k(randi(poolSize, nToTry, 1));
%     end

      objsToTry(n, X{2}(prevResponses)) = 0;

%     objsToTry(n, find(possibleObjs)) = 1;
end

if strcmp(method,'AS') || strcmp(method,'HFAS') || strcmp(method, 'HFQuick')

    %perform clustering
    nClust = round(size(P,1)/20); %maximum number of clusters we want to deal with
    display(['nClust=' num2str(nClust)]);
    if nClust > 40
        nClust = 40;
    end
    if nClust < nToAllocate
        nClust = nToAllocate;
    end

    chosenDocs = zeros(nK,nClust);
    
    diffAgents = 1:nK-1;
%      diffAgents = (1:nK)'; %use this if you want all the agents to get
%     the same objects to test against

    if exist(centroidsFile,'file')
        load(centroidsFile);
        display(['Last cluster update: ' num2str(clLastUpdate)]);
    else
        centroids = cell(1, size(diffAgents,2));
        clLastUpdate = 0;
    end
    % D = sparse(X{2}, X{1}, X{3}, size(Xfeat,2), nAllWorkers); %don't put this
    % in as it is reflected in P and has missing entries (gaps count as zeros
    % not as missing entries
    
    if length(X{1})-clLastUpdate>20 || clLastUpdate==0
        redoClusters = true;
        D = Xfeat';
        D = [D P];
    else
        redoClusters = false;
    end
        
    if ~iscell(centroids) || length(centroids)<size(diffAgents,2)
        centroids = cell(1, size(diffAgents,2));        
    end

    iter=1;
    for k=diffAgents
        
        centroids_i = centroids{iter};
        options = statset('UseParallel', 'always');       
        %redo clusters once every 20 iterations
        if redoClusters && ~isempty(centroids_i) && size(centroids_i,1)==nClust
            [clIdx, centroids_i, ~, ~] = kmeans(D, nClust, 'Start', centroids_i, 'EmptyAction', 'singleton', 'options', options);
            clLastUpdate = length(X{1});  
        elseif redoClusters
            [clIdx, centroids_i, ~, ~] = kmeans(D, nClust, 'EmptyAction', 'singleton', 'options', options);
            clLastUpdate = length(X{1});  
        end
        centroids{iter} = full(centroids_i);

        clustIds = unique(clIdx);
        nClust = length(clustIds);

        %do sub-clustering of each cluster so that each worker is tested against
        %different subsets -> better coverage of input space esp. when choosing
        %between workers
        for cl=1:nClust
            members = find(clIdx==clustIds(cl));
            memberIdx = randi(length(members),1,1); %select randomly from each cluster
            chosenDocs(k, cl) = members(memberIdx);
            memberIdx = 1;    
            while sum(objsToTry(k, chosenDocs(1, cl))==0)>0 && memberIdx<=length(members)
                chosenDocs(k, cl) = members(memberIdx);
                memberIdx = memberIdx+1;
            end
        end
    
        iter = iter+1;
    end
  
    save(centroidsFile, 'centroids', 'clLastUpdate', 'clIdx');
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

IG = zeros(nK, length(objsToTry)) - inf;

if strcmp(method,'AS') || strcmp(method,'HFAS') %complete IG table using all proposed documents from the clustering step
    
     if strcmp(method,'AS') || sum(chosenDocs(nK,:))==0
         nToCompare = nK-1;
     else
         nToCompare = nK;
     end
    
     for n=1:nToCompare
        display(['evaluating agent ' num2str(K(n))]);
        k = K(n);
        
        IGn = zeros(1, length(chosenDocs(n,:)));
        Pin = Pi(:,:,n);
        lnPin = lnPi(:,:,n);
        objs_n = objsToTry(n, chosenDocs(n,:));
        parfor i=1:length(chosenDocs(n,:))
            if objs_n(i)==0
                IGn(i) = -inf;
                continue
            else
                %display([num2str(i)]);
                if IGTargetsOnly
                    IGn(i) = evaluatePair1(k, chosenDocs(n,i), P, logJoint, Pin, lnPin, X, ETnow, combiner);
                else
                    IGn(i) = evaluatePair3(k, chosenDocs(n,i), P, logJoint, Alpha_all, AlphaFeat, Nu, Pin, lnPin, X, ETnow, combiner);
                end
                
            end
        end
        
        IG(n,chosenDocs(n,:)) = IGn;
     end   
%      if ~strcmp(method,'AS')
%          priorHt = -sum(P(chosenDocs(1,:),:) .* log(P(chosenDocs(1,:),:)),2);
%          [~, mostUncertain] = max(priorHt);    
%          MU = chosenDocs(1,mostUncertain);
%          display(['Unknown worker for most uncertain: ' num2str(MU)]);
%          if IGTargetsOnly
%             IG(nK,MU) = evaluatePair1(K(end), MU, P, logJoint, Pi(:,:,nK), lnPi(:,:,nK), X, ETnow, combiner);
%         else
%             IG(nK,MU) = evaluatePair3(K(end), MU, P, logJoint, Alpha_all, AlphaFeat, Nu, Pi(:,:,nK), lnPi(:,:,nK), X, ETnow, combiner);        
%         end
%      end
elseif strcmp(method, 'HFQuick') %complete IG table with only one doc as a tester

    P(chosenDocs(1,1),:)

    priorHt = -sum(P(chosenDocs(1,:),:) .* log(P(chosenDocs(1,:),:)),2);
    [~, mostUncertain] = max(priorHt);
    MU = chosenDocs(1,mostUncertain);

    for n=1:nK%(nK-1)
        display(['evaluating agent ' num2str(K(n))]);
        k = K(n);
        i = MU;
        if objsToTry(n, MU)==0
            i = find(objsToTry(n,:), 'first');
        end
%         display(num2str(i));
        if IGTargetsOnly
            IG(n, :) = evaluatePair1(k, i, P, logJoint, Pi(:,:,n), lnPi(:,:,n), X, ETnow, combiner);
        else 
            IG(n, :) = evaluatePair3(k, i, P, logJoint, Alpha_all, AlphaFeat, Nu, Pi(:,:,n), lnPi(:,:,n), X, ETnow, combiner);        
        end
    end
    
elseif strcmp(method, 'ASEntropy')
    %complete IG table using maximum entropy
    IG = repmat(-sum(P .* log(P),2)', size(IG,1), 1);   
elseif strcmp(method, 'Random2')
    IG = rand(size(IG));
    keepNewWorker = false;
    for n=1:nK
	k = K(n);
	if sum(X{1}==k)>20
	    IG(n,:)=0;
	    keepNewWorker = true;
	end
    end
    if ~keepNewWorker
	IG(nK,:) = 0; %don't allow the unkown worker to be assigned! No hiring and firing should take place.
    end
display('How do we choose the number after which we fire people? If we demonstrate differences between dynamic changes of individuals, (e.g. some improve, some get worse, some are good for ages) this is obviously a poor strategy');
elseif strcmp(method, 'Random')
    IG = rand(size(IG));
    IG(nK,:) = 0; %don't allow the unkown worker to be assigned! No hiring and firing should take place.
    
elseif strcmp(method, 'OS')
    %Do online screening: assign each worker a confidence level, fill in IG matrix with this.
    %Don't differentiate between docs though.
    minPrevResp = inf;
    display('trust:');
    for n=1:nK
        trust = singleTrustValue(Alpha(:,:,n));
        IG(n,:) = trust;

        prevResponses = sum(X{1}==K(n));
        if minPrevResp > prevResponses && n<nK
            minPrevResp = prevResponses;
        end
        
    end

    IG(:,1)
end

for n=1:nK
    IG(n, X{2}(X{1}==K(n))) = 0;
end

%select the optimal pairs of workers/documents.
%exclude those who are more useless than a random new member.
taskIdxs = zeros(1, nToAllocate);
permittedWorkers = zeros(1, nToAllocate);

remainIG = IG;

%ensure each worker has a task if better than the unknown worker
nExtraTasks = nToAllocate - (nK-1);
display(['n extra tasks ' num2str(nExtraTasks)]);

newbies = [];
for k=1:nK-1
    if sum(X{1}==K(k))<5
        newbies = [newbies k];
    end
    
    display(['remainIG for k=' num2str(k) ', number of candidate docs=' num2str(sum(remainIG(k,:)>0))]);
end

display(['Have to allocate to the newbies: ' num2str(newbies)]);


%find most certain object; assign to the newbies if we exclude someone
%  priorHt = -(P(chosenDocs(1,:),:) .* log(P(chosenDocs(1,:),:)),2);

% [kMax kMaxIdx] = max(remainIG, [], 2);
% display(['Max IG for each agent: ' num2str(kMax')]);
% [sortedKmax sortedKmaxIdx] = sort(kMax, 1, 'descend');

agentScores = sparse(1,nAllWorkers);
unknownScores = zeros(1, nAllWorkers);

for n=1:nToAllocate
    
    [kMax kMaxIdx] = max(remainIG, [], 2);
    
%     if exist('igSet.txt','file') && n==1
%         if exist('igfile.mat','file')
%             load('igfile.mat', '-ascii');
%         else
%             ig_set = zeros(1,100);
%         end
%         agentIds = K;
%         agentIds(nK) = 100;
%         ig_set = [ig_set; sparse(1,agentIds,kMax)];
% 
%         save('igfile.mat','ig_set', '-ascii');
%     end
    
    display(['Max IG for each agent: ' num2str(kMax')]);
    [overallMax overallMaxIdx] = max(kMax);
    
    agentScores(K(overallMaxIdx)) = overallMax;
    
%     overallMax = sortedKmax(n);
%     overallMaxIdx = sortedKmaxIdx(n);
    
    if overallMax==0 %unlikely case that we run out of tasks to allocate
        display('this allocation would give zero EIG');
        newIdx = overallMaxIdx;
        while ismember(newIdx, permittedWorkers) 
            newIdx = newIdx +1;
            if  newIdx>size(IG,2)
                newIdx = 1;
            elseif newIdx==overallMaxIdx
                display('we allocated absolutely everything, so just going to allocate documents multiple times at random to new workers');
            end
        end
        overallMaxIdx = newIdx;
    end
    
    while ~isempty(newbies) && ~ismember(overallMaxIdx, newbies) && nToAllocate-n < length(newbies)
        display(['getting the newbies in; will have to fire ' num2str(overallMaxIdx)]);
        remainIG(overallMaxIdx, :) = -inf;
        [kMax kMaxIdx] = max(remainIG, [], 2);
        display(['new Max IG for each agent: ' num2str(kMax')]);
        [overallMax overallMaxIdx] = max(kMax);
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
    
    k = K(overallMaxIdx);
    
    theChosenOne = overallMaxIdx;
    
    if (strcmp(method,'HFAS') || strcmp( method, 'HFQuick')) && ~ismember(overallMaxIdx,newbies)
        %Methods that allow the unknown worker to be assigned instead of
        %the proposed worker (firing). OS method would already do this when
        %finding the max, since unknown worker is calculated in advance
        if IG(nK,i)<=0%unknown worker has not been evaluated
            if strcmp( method, 'HFQuick') 
                i = MU;
            end
            if IGTargetsOnly
                IG(nK,i) = evaluatePair1(K(end), i, P, logJoint, Pi(:,:,nK), lnPi(:,:,nK), X, ETnow, combiner);
            else
                IG(nK,i) = evaluatePair3(K(end), i, P, logJoint, Alpha_all, AlphaFeat, Nu, Pi(:,:,nK), lnPi(:,:,nK), X, ETnow, combiner);    
            end
            remainIG(nK,i) = IG(nK,i);
        end
        
%         display(['unknown worker IG: ' num2str(remainIG(nK,i)) ', ' num2str(i) ' vs. ' num2str(remainIG(overallMaxIdx,i))]);

        if sum(X{1}==k) > 5 && remainIG(nK,i) - overallMax > 0.0001
            display(['diff: ' num2str(remainIG(nK,i) - overallMax)]);
            overallMax = remainIG(nK,i);
            overallMaxIdx = nK;
            %chosen the unknown worker - don't give tasks to the others
            k = 0;
        elseif overallMaxIdx>nAllWorkers
            k = 0;
        end        
        
        Alpha(:,:,overallMaxIdx)
    else %method = AS, ASEntropy, Random; i.e. methods that don't fire workers online
        k = K(overallMaxIdx);
    end
    
%     if k==0 || k==nAllWorkers+1
%         [pMax, silverIdx] = max(max(P(:, 2:end), [], 2), [], 1);
%         P(silverIdx,:) = 0;
%         remainIG(
%         i = silverIdx;
%         display(['Silver task for new workers: ' num2str(silverIdx) ' with p(most likely topic)=' num2str(pMax)]);
%     end

    unknownScores(K(theChosenOne)) = IG(nK,i);
        
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
            remainIG(overallMaxIdx,:) = -inf;
        end
        %remove from list of waiting newbies
        newbies(newbies==overallMaxIdx) = [];
    end
    
    permittedWorkers(n) = k;
    
    remainIG(:,i) = -inf; %remove the current object and look for new tasks
end

if size(K,1)>size(K,2)
    K = K';
end
for k=K
    if ~ismember(k, permittedWorkers) && k<=nAllWorkers
        display(['Excluding ' num2str(k)]);
        Alpha(:,:,K==k)
        excludedWorkers = unique([excludedWorkers; k]);
    end
    
%     n = find(K==k);
     
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
    
    trust = round(trust.*1000) ./ 1000;
end