function IG = batchInfoGain(i, state, calcPiInfoGain, simplify)
%IMMEDIATEINFOGAIN Expected information gain for assigning i to k. Doesn't care
%about future steps.
%Knockon is set to true if ommitted, meaning we also include the
%information gained about other classifiers and targets
    if nargin < 3
        calcPiInfoGain = false;
    end

    %the prior
    prET = state{3};
    prAlpha = state{1};
            
    nAgents = size(prAlpha,3);
    if nargin < 4
        simplify = true; %use the simpler, rougher calculation?
        %if using simpler calculation, shall we add on IG about pi as well
        %as IG about future data points? No, the value of IG about pi should
        %be taken into account by the knock-on effect on IG of other
        %labels.
    end
    pCgivT = prAlpha ./ repmat(sum(prAlpha,2), [1,size(prAlpha,2),1]);
% infgain.pClassifierOutput(state, i, 0, simplify);

    
    IG = 0;
    if ndims(prET)==2
        pT = repmat(prET(:,i), [1, size(pCgivT,2),nAgents]);%repmat(sum(pCT,2), 1, size(pCT,2)); %prior - aren't these two the same???
    else
        pT = repmat(prET(:,i,:), [1, size(pCgivT,2)]);
    end
    pCT = pCgivT .* pT;
    pC = repmat(sum(pCT,1), size(pCT,1), 1);
    pTgivC = pCT./pC;
    objectIG = sum(sum(pCT.*log(pTgivC) - (pC.*pT.*log(pT)), 2), 1);  
    %postAlpha is a cell array of alpha objects, each with 3
    %dimensions, one for each possible score
    for l=1:size(prAlpha,2)
        postAlpha = prAlpha;
        postAlpha(:,l,:) = postAlpha(:,l,:) + pTgivC(:,l,:);  

        if nargin > 2 && calcPiInfoGain
            postNu = repmat(state{2}', [1,1,nAgents]) + pTgivC(:,l,:) - pT(:,l,:);
            IG = IG + piInfoGain(prAlpha, postAlpha, state{2}, postNu);             
        end     
        
        % *** Also need to calculate knowck-on effect on other t?
        % Can we do a joint IG calculation? Still need to add on an extra
        % term for future benefit?
    end
    
    IG = IG + reshape(objectIG, [1,nAgents]);
end

function IG = piInfoGain(prAlpha, postAlpha, prNu, postNu)
%not the information gain about the confusion matrices. This is really about value of the confusion
%matrices.
    [hdRating hdClassRatings igRating] = infgain.piRating(prAlpha, prNu);
    [hdRating hdClassRatings postIgRating] = infgain.piRating(postAlpha, postNu);

    IG = postIgRating - igRating;
end