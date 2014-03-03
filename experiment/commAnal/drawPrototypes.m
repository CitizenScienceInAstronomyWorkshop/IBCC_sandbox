function drawPrototypes(g, P, Pi)
figure;

n=0;
[A, C] = max(P, [], 1);
[membVals discreteComms] = max(P, [], 2);
for c=1:size(P,2)
    
    commAgents = find(discreteComms==c);
    if numel(commAgents) < 1
        display(['skipping ' num2str(c)]);
        continue;
    else
        n = n + 1;
    end
    
    agent = C(c);
    
    PiComm = Pi(:, :, agent);
    
    avgPi = PiComm;
    
    subplot(1, ceil(numel(g)), n);
    
    X = [0 1];
    bar3(X, avgPi);
end
end