nPoints = 10001;

X = 0:nPoints-1;
X = X ./ (nPoints-1);

A = log(X) - log(1-X);
                
C = ones(1, nPoints) .* 1;
M = X;
contributions = normpdf(ones(nPoints, 1)*X,...    
                    M' * ones(1,nPoints),...
                    C' * ones(1, nPoints));
                
                
fX = sum(contributions, 1) ./ nPoints;
                
figure;
plot(A, fX);

hold all