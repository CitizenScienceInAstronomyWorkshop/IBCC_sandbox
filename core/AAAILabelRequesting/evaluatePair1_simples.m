function [IG Hnow Hnew] = evaluatePair1(k, i, P, logJoint, Pi, lnPi, X, ETnow, combiner)
    % Practically, we can skip pairs where:
    %   1. Object and worker are both well known, unless there are lots
    %   of similar objects or lots of objects classified by worker. 
    %   2. Object and worker are both poorly known unless there are
    %   lots of similar objects.
    nScores = size(Pi,2);
    pCi = sum(Pi .* repmat(P(i,:)',1,nScores), 1);
 
    %simulate the updates                 
%     logJoint(isinf(logJoint)) = 0;
%     blurgh = sum(repmat(pCi,size(Pi,1),1) .* lnPi, 2);
%     logJoint(:,i) = logJoint(:,i) + blurgh;
    
    Hnow = -sum(P' .* log(P'),1);%-sum(P' .* (logJoint), 1) + log(sum(exp(logJoint),1));
    
    Hnew = zeros(nScores, size(P,1));

    for c=1:nScores
        
        X{1} = [X{1}; k];
        X{2} = [X{2}; i];
        X{3} = [X{3}; c];
%         combCopy = combiners.bcc.MixedIbccVb(bccSettings, nFeat, featSettings, ...
%             Xfeat, nAllWorkers, 1, Tvec, [], nClasses, nClasses);
%         combCopy.initAll = true;
        combiner.combinedPost = ETnow;
        Pc = combiner.combineDecisions(X);   
        
%         lnPtc = combiner.logJoint;
%         lnPtc(isinf(lnPtc)) = 0;
        
        Hnew(c,:) = -sum(Pc .* log(Pc), 1); %-sum(Pc .* lnPtc, 1) +log(sum(exp(lnPtc),1));
        Hnew(c,:) = Hnew(c,:) .* pCi(c);

        X{1} = X{1}(1:(end-1));
        X{2} = X{2}(1:(end-1));
        X{3} = X{3}(1:(end-1));
    end
%     Hnow
%     Hnew
%     pCi
%     Hnew .* pCi
    
    Hnew = sum(Hnew, 1);
%     Hnew

    IG = Hnow - Hnew;
    
%     display(['Number of negative IG points: ' num2str(sum(IG<0)) ' and their value: ' num2str(sum(IG(IG<0))) '; as a fraction of total IG: ' num2str(sum(IG(IG<0))./sum(IG))]);
    
    IG = sum(IG);
end

