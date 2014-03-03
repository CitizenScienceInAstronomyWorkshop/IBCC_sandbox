function [clusterMembershipChanges] = plotResultsClusterColours( a, ...
    agentClusterMembership, results, colourSet, nSamples, ...
    comparison_line )
% Plots the results for an agent onto the current figure and subplot,
% colouring the line by cluster membership, and drawing points where
% changes in cluster membership occur.
%Colour_set should be a set of colours corresponding to clusters.
%Comparison_line is a boolean flag indicating whether this should be a plot
%for comparison, in which case it is widened and dashed to distinguish from 
% other lines in the graph.


    clusterMembershipChanges = {};
    previous_cluster = 0;
    for n=1:nSamples
        if previous_cluster ~= agentClusterMembership(a, n);
            clusterMembershipChanges{ ...
                length(clusterMembershipChanges)+1 } = n;

            previous_cluster = agentClusterMembership(a, n);
        end

    end

    start_idx = clusterMembershipChanges{1};

    for i=2:length(clusterMembershipChanges)+1
        if i > length(clusterMembershipChanges)
            c = nSamples;
        else
            c = clusterMembershipChanges{i};
        end

        colour = colourSet(agentClusterMembership(a, start_idx), :);

        plot(start_idx, results(a, start_idx), 'o', 'Color', colour, 'MarkerSize', 4);

        if exist('comparison_line', 'var') && comparison_line
            %lighten the colour
            colour = colour + [0.6 0.6 0.6];
            colour(colour<0) = 0;
            colour(colour>1) = 1;
            plot(start_idx:c, results(a, start_idx:c), '--', 'Color', colour, 'LineWidth', 2);
        else
            plot(start_idx:c, results(a, start_idx:c), 'Color', colour, 'LineWidth', 1); 
        end

        start_idx = c;

        hold on

    end
end

