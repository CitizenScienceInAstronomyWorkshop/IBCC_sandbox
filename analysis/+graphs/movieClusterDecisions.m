%run cluster_agents_kmeans first

colourSet = createColourSet(K);

figure;
for n=1:nSamples
    hold off
    
    for k=1:K
        cluster_agents = find(agentClusterMembership(:, n)==k);
        cluster_results = baseLogodds(agentClusterMembership(:, n)==k, n);
        
        plot(cluster_results, cluster_agents, 'x', ...
            'MarkerSize',10, 'Color',colourSet(k, :));
        
        title(sprintf('n=%i, changepoints: %i -> inf. %i -> useless', ...
            n, ...
            labelledTestData(n, nSensors+2), ...
            labelledTestData(n, nSensors+3) ...
            ));
        
        hold all
    end
    %M(n) = getframe;
    %drawnow;
    pause;
end
%movie(M, 1, 1);