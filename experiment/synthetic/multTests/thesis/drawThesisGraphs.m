function [aucs, meanAucs, devAucs, agentAucs] = drawThesisGraphs( nDatasets, expSet, settingVals, varName,...
    combinedPostKnowns, agentPostKnowns, labelsKnowns, combMethods, nSettings, ...
    graphColors, graphMarkers, graph1, graph2,...
    graph1Colors, graph1Markers, graph2Colors, graph2Markers, specialExp)

    aucs = zeros(nSettings, length(combMethods)+2, nDatasets);
    totalAgents = zeros(nSettings,1);
    figs = zeros(nSettings+1, 1);
    
    agentAucs = {};
    
    for i=1:nSettings
        testLabels = [];
        avgAgentAucs = 0;
        combinedPost = [];
        
        agentAucsI = [];
        
        for d=1:nDatasets
            r=1;
            combinedPost = [combinedPost combinedPostKnowns{d,i}(r,:,:)];
            testLabels = [testLabels labelsKnowns{d,i}(r,:,:)];
            
            aucs(i,1:end-2, d) = aucs(i,1:end-2, d) + graphs.ClassifierPerformanceGraph.drawRoc(...
                               reshape(combinedPostKnowns{d,i}(r,:,:),size(combinedPostKnowns{d,i},2),length(combMethods))',...
                               labelsKnowns{d,i}(r,:,:)-1, combMethods,false,true);

            agentAucs_d = graphs.ClassifierPerformanceGraph.drawRoc(...
                            agentPostKnowns{d,i}, labelsKnowns{d,i}(r,:,:)-1, cell(1,size(agentPostKnowns{d,i},1)), false, true);
            
            aucs(i,end, d) = max(agentAucs_d);
            
            agentAucsI = [agentAucsI agentAucs_d];
            
            aucs(i,end-1, d) = sum(agentAucs_d)./length(agentAucs_d);
            
%              aucs(i,end) = aucs(i,end) + sum(agentAucs);
%              totalAgents(i) = totalAgents(i) + length(agentAucs);            
        end
        
        agentAucs{length(agentAucs)+1} = agentAucsI;
        
        testLabels = testLabels -1;
        combinedPost = reshape(combinedPost, size(combinedPost,2), length(combMethods))';

        separateMethods = false;
        aucOnly = false;
        newFigure = false; %create our own handle
%         if i==nSettings %legend only on last graph
%             noLegend = false;
%         else
%             noLegend = true;
%         end
        noLegend = false;

        if ~isempty(graph2)
            figs(i) = figure('Position', [200, 400, 825, 300], 'PaperType', 'A4', 'PaperPositionMode', 'auto');
            subplot(1,2,1);
        else
            figs(i) = figure('Position', [200, 400, 450, 300], 'PaperType', 'A4', 'PaperPositionMode', 'auto');
        end
        
%          [aucs(i,graph1), ~, ~, lineHandles] 
        [~,~,~,lineHandles] = graphs.ClassifierPerformanceGraph.drawRoc(...
            combinedPost(graph1,:), testLabels, combMethods(graph1), separateMethods, aucOnly, newFigure, noLegend);        
        title(['ROC with ' varName '=' num2str(settingVals(i))]);
        hold all
        plot([0 1], [0 1], 'Color', 'black', 'linewidth', 1);
        grid on
        
        for i=1:length(lineHandles)
            %don't plot all the markers as there are too many.       
%             if ~strcmp(graph1Markers{i}(end), '-')
%                 set(lineHandles(i), 'Marker', graph1Markers{i}(end));
%             end
            set(lineHandles(i), 'Color', graph1Colors{i});
        end
        
        if exist('specialExp','var') && ~isempty(specialExp)
            set(lineHandles(specialExp), 'LineStyle', '-.');
            set(lineHandles(specialExp), 'Marker', '+');
        end
        
%         set(lineHandles(end), 'LineStyle', '-.');
%         set(lineHandles(end), 'Marker', 'x'); 
%         set(lineHandles(end), 'Color', 'magenta');
        
        if ~isempty(graph2)
            subplot(1,2,2);
%              [aucs(i,graph2), ~, ~, lineHandles] 
            [~,~,~,lineHandles] = graphs.ClassifierPerformanceGraph.drawRoc(...
                combinedPost(graph2,:), testLabels, combMethods(graph2), separateMethods, aucOnly, newFigure, noLegend);        
            %title(['ROC with ' varName '=' num2str(settingVals(i))]);
            hold all
            grid on
            plot([0 1], [0 1], 'Color', 'black', 'linewidth', 1);
            for i=1:length(lineHandles)
                %don't plot all the markers as there are too many.
%                 if ~strcmp(graph2Markers{i}(end), '-')
%                     set(lineHandles(i), 'Marker', graph2Markers{i}(end));
%                 end
                set(lineHandles(i), 'Color', graph2Colors{i});
            end            
        end 
    end
    
    meanAucs = sum(aucs,3) ./ nDatasets;
    varAucs = sum((aucs-repmat(meanAucs, [1,1,nDatasets])).^2, 3) ./ nDatasets;
    devAucs = varAucs.^0.5;
    
    combMethods{length(combMethods)+1} = 'Mean Individual';    
    combMethods{length(combMethods)+1} = 'Best Individual';    

    if ~isempty(graph2)
        figs(nSettings+1) = figure('Position', [200, 400, 1155, 420], 'PaperOrientation', 'landscape', 'PaperType', 'A4', 'PaperPositionMode', 'auto');
        subplot(1,2,1);
    else
        figs(nSettings+1) = figure('Position', [200, 400, 630, 420], 'PaperOrientation', 'landscape', 'PaperType', 'A4', 'PaperPositionMode', 'auto');
    end
    
    lineHandles = [];
    i=0;
    for yData=meanAucs(:,[graph1 end])
        i=i+1;
        lineHandles = [lineHandles plot(settingVals, yData, graph1Markers{i}, 'linewidth', 2,  'Color', graph1Colors{i})];
        hold all
    end
    if exist('specialExp','var') && ~isempty(specialExp)
        set(lineHandles(specialExp), 'LineStyle', '.');
    end
%      set(lineHandles(end-1), 'LineStyle', '-.');
%      set(lineHandles(end-1), 'Marker', 'x'); 
%      set(lineHandles(end-1), 'Color', 'black');
    
    set(lineHandles(end), 'LineStyle', '--');
    set(lineHandles(end), 'Marker', '+'); 
    set(lineHandles(end), 'Color', 'black');
    
    title(['AUC of Combined Decision Against ' varName]);
    xlabel(varName);
    ylabel('AUC');    
    legend(combMethods([graph1 end]), 'Location', 'Best');
    grid on
    if length(settingVals)>1
        axis([settingVals(1), settingVals(end), 0.5, 1]);
    end
    
    if ~isempty(graph2)
        subplot(1,2,2);
        lineHandles = [];
        i=1;
        for yData=meanAucs(:,[graph2 end])
            lineHandles = [lineHandles plot(settingVals, yData, graph2Markers{i}, 'linewidth', 2, 'Color', graph2Colors{i})];
            hold all
            i=i+1;
        end
        
%          set(lineHandles(end-1), 'LineStyle', '-.');
%          set(lineHandles(end-1), 'Marker', 'x'); 
%          set(lineHandles(end-1), 'Color', 'black');
    
        set(lineHandles(end), 'LineStyle', '--');
        set(lineHandles(end), 'Marker', '+'); 
        set(lineHandles(end), 'Color', 'black');        
        
        hold all;
%          if size(settingVals,2)<size(settingVals,1)
%              settingVals = settingVals';
%          end
%          plot(repmat(settingVals, 2, 1), meanAucs(:, [graph2 end])-devAucs(:,[graph2 end])
        
        %title(['AUC of Combined Decision against ' varName]);
        xlabel(varName);
        ylabel('AUC');
        legend(combMethods([graph2 end]), 'Location', 'Best');
        grid on
        axis([settingVals(1), settingVals(end), 0.5,  1]);
    end
    %save the figures
    for i=1:nSettings+1
        saveas(figs(i), [expSet.outputDir '/' num2str(i) '.pdf']);
        saveas(figs(i), [expSet.outputDir '/' num2str(i) '.fig']);
    end
end

