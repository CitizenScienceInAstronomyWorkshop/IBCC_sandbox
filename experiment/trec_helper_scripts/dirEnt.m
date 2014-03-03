function H = dirEnt( Alpha )
%DIRENT Calculate entropy of a dirichlet with parameter vector Alpha
B = prod(gamma(Alpha)) ./ gamma(sum(Alpha));
K = length(Alpha);
A0 = sum(Alpha);
H = - log(B) - (A0-K)*psi(A0) + sum( (Alpha-1).*psi(Alpha) );
end

