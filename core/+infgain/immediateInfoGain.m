function IG = immediateInfoGain(i, k, state, knockon, bcc, calcPiInfoGain)
%IMMEDIATEINFOGAIN Expected information gain for assigning i to k. Doesn't care
%about future steps.
%Knockon is set to true if ommitted, meaning we also include the
%information gained about other classifiers and targets
    if nargin < 4
        knockon = false;
        bcc = [];
        calcPiInfoGain = false;
    end
    IG = objectInfoGain(i,k,state, knockon, bcc, calcPiInfoGain);
end

function IG = objectInfoGain(i,k,state, knockon, bcc, calcPiInfoGain)
%KL divergence between the prior belief in class of i and the posterior given the response
%from k.  

    % p(t_i | c^k) = kappa_i*gamma(alpha+1)gamma(sum(alpha))/gamma(sum(alpha)+1)gamma(alpha)
    % D_kl(post||prior) = sum_j{ p(t_i|c^k) log{p(p(t_i|c^k) / p(t_i)}}
    %                   = sum_j{ p(t_i|c^k) + gammaln(alpha+1) + gammaln
    %                      (sum(alpha)) - gammaln (sum(alpha)+1) -
    %                      gammaln(alpha)
    
    %TASKS: 
    % 1. update iterations using vb
    % 2. info gain for the confusion matrices
    % 3. tidy up: don't need to separate the functions?
    % 4. add a switch or separate the code so that piInfoGain can be
    % avoided (part of 3?).
    
    %the prior
    prET = state{3};
    prAlpha = state{1};
            
    pCT = pClassifierOutput(state, i, k);
    pC = repmat(sum(pCT,1), size(pCT,1), 1);
    
    IG = 0;
    
    %recursive part
    if knockon
        bcc.C{1} = [bcc.C{1} k];
        bcc.C{2} = [bcc.C{2} i];
        origScores = bcc.C{3};
        
        pTgivC = zeros(size(prET,1), size(pCT,2), size(prET,2));
        pT = repmat( reshape(prET, [size(prET,1), 1, size(prET,2)]), 1, size(pCT,2));%prior
        
        for l=1:size(pCT,2)
            bcc.C{3} = [origScores l];
            %need to update pTgivC and postAlpha here - other values already calculated
            %postET = IBCCVB...
            [lnPi, lnP, postET, nIt, postAlpha, postNu] = bcc.iterate(nItPrev, lnPi, lnP, C, nAssets, ET, PostAlpha, agentIdx);
            pTgivC(:, l, :) = reshape(postET, [size(postET,1), 1, size(postET,2)]);
            if nargin > 5 && calcPiInfoGain
                IG = IG + piInfoGain(prAlpha, postAlpha, state{2}, postNu);
            end              
        end
        nUpdated = size(prET,2);
    else    
        pT = repmat(prET(:,i), 1, size(pCT,2));%repmat(sum(pCT,2), 1, size(pCT,2)); %prior - aren't these two the same???
            
        pTgivC = pCT./pC;
        nUpdated = 1;
        %postAlpha is a cell array of alpha objects, each with 3
        %dimensions, one for each possible score
        for l=1:size(prAlpha,2)
            postAlpha = prAlpha;
            postAlpha(:,l,k) = postAlpha(:,l,k) + pTgivC(:,l);  
            
            postNu = state{2} + pTgivC(:,l)' - pT(:,l)';
            
            if nargin > 5 && calcPiInfoGain
                IG = IG + piInfoGain(prAlpha, postAlpha, state{2}, postNu);             
            end            
        end
    end
    objectIG = sum(sum(pC .* pTgivC .* log(pTgivC ./ pT), 2), 1);
    if nUpdated>1
        objectIG = sum(objectIG,3);
    end    
    IG = IG + objectIG;
        

end

function IG = piInfoGain(prAlpha, postAlpha, prNu, postNu)
%not the information gain about the confusion matrices. This is really about value of the confusion
%matrices.
    [hdRating hdClassRatings igRating] = piRating(prAlpha, prNu);
    [hdRating hdClassRatings postIgRating] = piRating(postAlpha, postNu);

    IG = sum(postIgRating - igRating);
end