function [IG Hnow Hnew] = evaluatePair3(k, i, P, logJoint, Alpha, AlphaFeat, Nu, Pi, lnPi, X, ETnow, combiner)
    %like evaluatePair1 except we look at IG over all model variables, not just the targets.
    %This allows us to re-explore areas of the target space that would inform the confusion matrix.
    
    nScores = size(Pi,2);
    pCi = sum(Pi .* repmat(P(i,:)',1,nScores), 1);
 
    %simulate the updates                 
    logJoint(isinf(logJoint)) = 0;
    
    blurgh = sum(repmat(pCi,size(Pi,1),1) .* lnPi, 2);
    
    logJoint(:,i) = logJoint(:,i) + blurgh;
    
    %target label entropy
%      display('Hnow_t');
    Hnow_t = -sum(sum(P' .* (logJoint), 1), 2) + sum(log(sum(exp(logJoint),1)), 2);
    %class proportions entropy
%      display([ 'Hnow_kappa' num2str(Nu)]);
    Hnow_kappa = dirEntropy(Nu);
    %conf matrix entropy
%      display('Hnow_pi');
    Hnow_pi = sum(sum(dirEntropy(Alpha)));
    Hnow_piFeat = sum(sum(dirEntropy(AlphaFeat)));
    %Total current entropy
    %      display('Hnow');
    Hnow = Hnow_t + Hnow_kappa + Hnow_pi + Hnow_piFeat;
    
    Hnew = zeros(1,nScores);

    for c=1:nScores
        
        X{1} = [X{1}; k];
        X{2} = [X{2}; i];
        X{3} = [X{3}; c];

        combiner.combinedPost = ETnow;
        Pc = combiner.combineDecisions(X);   
        
        lnPtc = combiner.logJoint;
        lnPtc(isinf(lnPtc)) = 0;
        
        Hnew_ct =  - sum(sum( Pc .* lnPtc, 1), 2) + sum(log(sum(exp(lnPtc),1)), 2);
        %class proportions entropy
        Hnew_ckappa = dirEntropy(combiner.Nu);
        %conf matrix entropy
        Alpha = combiner.Alpha;
        if iscell(Alpha)
            Alpha = Alpha{1};
        end
        AlphFeat = combiner.featAlpha;
        Hnew_cpi = sum(sum(dirEntropy(Alpha))); 
        Hnew_cpiFeat = sum(sum(dirEntropy(AlphaFeat))); 

        Hnew(c) = Hnew_ct + Hnew_ckappa + Hnew_cpi + Hnew_cpiFeat;
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

