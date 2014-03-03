function [IG Hnow Hnew] = evaluatePair1(k, i, P, logJoint, Pi, X, ETnow, combiner)
    % Practically, we can skip pairs where:
    %   1. Object and worker are both well known, unless there are lots
    %   of similar objects or lots of objects classified by worker. 
    %   2. Object and worker are both poorly known unless there are
    %   lots of similar objects.
    nScores = size(Pi,2);
 
    %simulate the updates                 
    logJoint(isinf(logJoint)) = 0;
    Hnow = -sum(sum(P' .* logJoint, 1), 2) + sum(log(sum(exp(logJoint),1)), 2);

    pCi = sum(Pi .* repmat(P(i,:)',1,nScores), 1);

    Hnew = zeros(1,nScores);

    for c=1:nScores
        
        X{1} = [X{1}; k];
        X{2} = [X{2}; i];
        X{3} = [X{3}; c];
%         combCopy = combiners.bcc.MixedIbccVb(bccSettings, nFeat, featSettings, ...
%             Xfeat, nAllWorkers, 1, Tvec, [], nClasses, nClasses);
%         combCopy.initAll = true;
        combiner.combinedPost = ETnow;
        Pc = combiner.combineDecisions(X);   
        
        lnPtc = combiner.logJoint;
        lnPtc(isinf(lnPtc)) = 0;
        
        Hnew(c) =  - sum(sum( Pc .* lnPtc, 1), 2) + sum(log(sum(exp(lnPtc),1)), 2);
        
        Hnew(c) = Hnew(c);

        X{1} = X{1}(1:(end-1));
        X{2} = X{2}(1:(end-1));
        X{3} = X{3}(1:(end-1));
    end
%     Hnow
%     Hnew
%     pCi
%     Hnew .* pCi
    
    Hnew = sum(Hnew .* pCi);
    
%     Hnew

    IG = Hnow - Hnew;
end

