function plotBasic(agentIdx, results, colourSet, nSamples, width, comparison_line)

    colour = colourSet(agentIdx, :);

    if exist('comparison_line','var') && comparison_line
        %lighten the colour
        colour = colour + [0.6 0.6 0.6];
        colour(colour<0) = 0;
        colour(colour>1) = 1;
        
        %changed line style from '--
        plot(1:nSamples, results(agentIdx, 1:nSamples), '--', 'Color', colour, 'LineWidth', width+1);
    else
        %changed line style from default (not specified)
        plot(1:nSamples, results(agentIdx, 1:nSamples), 'Color', colour, 'LineWidth', width);
    end

    hold on    
end

