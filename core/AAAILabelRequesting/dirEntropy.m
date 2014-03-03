function H = dirEntropy(A)
    %entropy of a dirichlet-distributed variable with parameter vector A
    
    As = sum(A, 2);
    
    %beta of A
    lnBetaA = sum(gammaln(A), 2) - gammaln(As);
    
    K = size(A,2);
    
    H = lnBetaA + ((As-K) .* psi(As)) - sum( (A-1).*psi(A), 2 );
end