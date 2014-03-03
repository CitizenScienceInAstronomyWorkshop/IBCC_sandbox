distance_matrix;

hold off;
figure;
for n=1:nSamples
    
    for a=1:nAgents
        for b=1:nAgents
            plot(a, b, 'o','MarkerFaceColor', [1 S_agents(a, b, n)  0],'MarkerSize', S_agents(a, b, n)*20, 'MarkerEdgeColor', [1 1 0]);
            text(a, b, sprintf('%3.3f', S_agents(a, b, n) ) );
            xlabel('agent number');
            ylabel('agent number');
            axis([0 nAgents 0 nAgents]);
            hold on
        end
    end    
    M(n) = getframe;
    hold off
end

movie(M, 1, 1);