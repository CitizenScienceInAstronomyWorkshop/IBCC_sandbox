function pCT = pClassifierOutput(state, i, k, simplify)
%pClassifierOutput Probaility of the classifier outputs of k given the object i
% Returns a matrix of joint classifier output/class label probabilities

if nargin<3 || k<1
    nAgents = size(state{1}, 3);
    k = 1:nAgents;
else
    nAgents = 1;
end

if nargin < 4
    simplify = true;
end

Kappa = state{3}(:,i);
Alpha = state{1};

if simplify
    pCgivT = Alpha ./ repmat(sum(Alpha,2), [1,size(Alpha,2),1]);

    
    
else
    gammasumi = repmat(gammaln(sum(Alpha(:,:,k),2)), [1, size(Alpha,2), 1]);
    gammasumin = repmat(gammaln(sum(Alpha(:,:,k),2)+1), [1, size(Alpha,2), 1]);
    logpCgivT = gammasumi - gammasumin  + gammaln(Alpha(:,:,k)+1) - gammaln(Alpha(:,:,k));
    pCT = exp(logpCgivT) .* repmat(Kappa, [1, size(Alpha,2), nAgents]);
end

