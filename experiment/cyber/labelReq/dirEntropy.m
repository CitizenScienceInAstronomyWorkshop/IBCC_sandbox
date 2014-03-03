function h = dirEntropy(A)
lnB = sum(gammaln(A), 2) - gammaln(sum(A, 2));
Asum = sum(A,2);
K = size(A,2); %nScores
h = lnB + (Asum-K).*psi(Asum) - sum(A-1.*psi(A), 2);
end