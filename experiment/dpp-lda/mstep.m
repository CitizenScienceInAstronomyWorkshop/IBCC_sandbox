function [ betastar ] = mstep( phi, w, beta, lambda, a )

% a = hyperparameter for beta
% beta_kv = p( w_mn = v | z_mn = k )
% phi_mnk = p( z_mn = k | w_mn )
% w_mn^(v) = { 1 if w_mn = v, 0 if w_mn ~= v
% betastar_kv  = sum_M sum_N_m phi_mnk w_mn^(v) = count of word v occuring
% in topic k given a distribution over topics for each word occurrence /
% count of all words in topic k = p( w_mn = v | phi_mnk, w_mn )

K = size(beta, 1);
rho = 0.5;

%kernel prior calculation
A = repmat(a, K, 1);
pi_beta = prod( beta.^A, 2) .* gamma(sum(A,2)) ./ prod(gamma(A), 2);

R_beta = (beta.^rho)*(beta.^rho)' ./ (sqrt(sum(beta.^(2*rho), 2)) * sqrt(sum(beta'.^(2*rho), 1)));

k_beta = R_beta .* sqrt(pi_beta *pi_beta');
lndet_k_beta = log(det(R_beta)) + sum(pi_beta);
beta_prior = lambda .* k_beta;

%Find gradient of this with respect to beta. Subtract the gradient from
%beta.
% = trace( K_beta ^ -1 d K_beta/d beta



end
