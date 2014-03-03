
function avgPis = drawWeightedMeans(g, P, Pi)
figure;
n=0;
[membVals discreteComms] = max(P, [], 2);

avgPis = zeros(size(Pi,1), size(Pi,2), length(g));

for c=1:size(P,2)
    commAgents = find(discreteComms==c);
    if numel(commAgents) < 1
        continue;
    else
        n = n + 1;
    end
    
    PiComm = Pi(:, :, commAgents);
    for a=1:numel(commAgents)
        PiComm(:,:,a) = PiComm(:,:,a) .* P(a,c);
    end
    
    avgPi = sum(PiComm, 3);
    %normalise
    avgPi = avgPi ./ repmat(sum(avgPi, 2), 1, size(avgPi,2));
    avgPis(:,:,n) = avgPi;
    subplot(1, ceil(numel(g)), n);
    
    if size(avgPi, 1) == 2
        X = [0 1];
        bar3(X, avgPi);
    else
        bar3(avgPi);
    end
end
end