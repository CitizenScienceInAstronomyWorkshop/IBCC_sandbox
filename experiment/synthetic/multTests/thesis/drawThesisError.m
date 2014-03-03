function [results] = drawThesisError( yLabel, expSet, settingVals, varName, ...
    combinedError, agentError, meanAgentError, combMethods,graph1,graph2, ...
    graph1Colors, graph1Markers, graph2Colors, graph2Markers)

    results = [combinedError; meanAgentError; agentError];

    if ~isempty(graph2)
        fig = figure('Position', [300, 400, 1155, 420], 'PaperType', 'A4','PaperOrientation', 'landscape', 'PaperPositionMode', 'auto');
        subplot(1,2,1);
    else
        fig = figure('Position', [300, 400, 630, 420], 'PaperType', 'A4', 'PaperOrientation', 'landscape', 'PaperPositionMode', 'auto');
    end
    
    lineHandles = [];
    i=1;
    for yData=results([graph1 end],:)'
        lineHandles = [lineHandles plot(settingVals, yData, graph1Markers{i}, 'linewidth', 2, 'Color', graph1Colors{i})];
        hold all
        i=i+1;
    end
    
%      set(lineHandles(end-1), 'LineStyle', '-.');
%      set(lineHandles(end-1), 'Marker', 'x');       
%      set(lineHandles(end-1), 'Color', 'black');    
    
    set(lineHandles(end), 'LineStyle', '--');
    set(lineHandles(end), 'Marker', '+');       
    set(lineHandles(end), 'Color', 'black');
    
    
    if exist('specialExp','var') && ~isempty(specialExp)
        set(lineHandles(specialExp), 'LineStyle', '.');
    end
    grid on
    if length(settingVals)>1
        axis([settingVals(1), settingVals(end), 0, 0.7]);     
    end
    title([yLabel ' for Different Values of ' varName]);
    xlabel(varName);
    ylabel(yLabel);    
    
    
    combMethods{length(combMethods)+1} = 'Mean Individual'; 
    combMethods{length(combMethods)+1} = 'Best Individual';   
    legend(combMethods([graph1 end]), 'Location', 'Best');
    
    if ~isempty(graph2)
        subplot(1,2,2);
        lineHandles = [];
        i=1;
        for yData=results([graph2 end],:)'
            lineHandles = [lineHandles plot(settingVals, yData, graph2Markers{i}, 'linewidth', 2, 'Color', graph2Colors{i})];
            hold all
            i=i+1;
        end
            set(lineHandles(end), 'LineStyle', '--');
        set(lineHandles(end), 'Marker', '+');     
        set(lineHandles(end), 'Color', 'black');
        
%          set(lineHandles(end-1), 'LineStyle', '-.');
%          set(lineHandles(end-1), 'Marker', 'x'); 
%          set(lineHandles(end-1), 'Color', 'black');    
        grid on
        %title(['Bayes Error for Different Values of ' varName]);
        xlabel(varName);
        ylabel(yLabel);
        legend(combMethods([graph2 end]), 'Location', 'Best');
        axis([settingVals(1), settingVals(end), 0, 0.7]);        
    end

    %save the figures
    saveas(fig, [expSet.outputDir '/errors.pdf']);
    saveas(fig, [expSet.outputDir '/errors.fig']);
end

