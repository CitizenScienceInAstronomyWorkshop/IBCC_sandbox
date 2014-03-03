classdef ClassifierPerformanceGraph < handle
    %draw graphs for the output of a set of classifiers. A
    %ClassifierPerformanceGraph object contains settings for plotting
    %different graphs for the same set of classifiers. 
    
    properties
        expSettings
        
        %agent objects required if we want to plot informedness of agents
        agents
        
        %index of an agent/classifier to be used purely as
        %comparison/benchmark for the other agents
        compareIdx
        %flag - whether to use compareIdx or not
        compare
        
        %these may be empty if the classifiers have not been clustered.
        clusterMembership
        clusterChanges
        
        scrsz
        
        colourSet
        clusterStrings
        agentStrings
        
        %should the graphs show which clusters the classifiers are in?
        clustered
        
        %sizes, numbers, dimensions etc.
        A
        K
        N
        S
        Si
        
        %sensor changepoints
        changes
        plotChangesSimple = true;
    end
    
    methods (Static)
        function Errors = calculateBayesErrors(resultsPost, nSamples, correctLabels)
            correctLabels = correctLabels -1;
            if size(correctLabels, 2) > size(correctLabels, 1)
                correctLabels = correctLabels';
            end            
            
            nAgents = size(resultsPost, 1);
            classTable = correctLabels(1:nSamples, 1) * ones(1, nAgents);
                       
            Errors = abs( classTable' - resultsPost(:, 1:nSamples) );
        end
        
        % ROC CURVE
        function [auc,fig,thresholds] = drawRoc(labelledCombinedPost, knownLabels, combMethods, separateMethods, aucOnly, newFigure)
            if nargin < 4
                separateMethods = false;
            end
            nCombMethods = length(combMethods);
            
            if size(knownLabels,1) > 1
                knownLabels = knownLabels';
            end            
            
            if max(knownLabels) > 1
                nClasses = max(knownLabels) + 1; %is this plus one needed?
            else
                nClasses = 1;
            end            
           
            if nClasses > 2
                for j=1:nClasses
                    labelledCombinedPost{j} = labelledCombinedPost{j}(:, knownLabels>=0);
                end
            else
                labelledCombinedPost = labelledCombinedPost(:, knownLabels>=0);
            end
            knownLabels = knownLabels(knownLabels>=0);           
             
            auc = zeros(nClasses, length(combMethods));
            
            for j=0:nClasses-1
                
                if nClasses > 1
                    labels = knownLabels==j;
                    results = labelledCombinedPost{j+1};
                else
                    labels = knownLabels;
                    results = labelledCombinedPost;
                end
                
                [tpr fpr thresholds] = roc(repmat(labels, nCombMethods, 1), results);
                if nargin < 5 
                    fig=figure('Visible', 'Off');
                elseif ~aucOnly && newFigure
                    fig=figure;
                end
                for c=1:length(combMethods)

                    if separateMethods
                        subplot(2, ceil(length(combMethods)/2), c);
                    end

                    if iscell(tpr)
                        tpr_c = [tpr{c} 1];
                        fpr_c = [fpr{c} 1];
                    else
                        tpr_c = [tpr 1];
                        fpr_c = [fpr 1];
                    end
                    
                    auc_jc =trapz(fpr_c, tpr_c);
                    auc(j+1, c) = auc_jc;
                    
                    if nargin < 5 || ~aucOnly
                        plot(fpr_c, tpr_c, 'linewidth', 2, 'clipping', 'off');   
                        hold all;
                    
                        if separateMethods
                            legend(combMethods{c});
                        end
                    end
                end
                if nargin < 5 || ~aucOnly
                    xlabel('False Positive Rate');
                    ylabel('True Positive Rate');
                    
                    if ~separateMethods
                        legend(combMethods);            
                    end
                end
            end
        end        
    end
    
    methods
        function obj = ClassifierPerformanceGraph(expSettings, ...
                clusterMembership, compareIdx, legendStrings, changes, agents)
            
            if exist('agents','var')
                obj.agents = agents;
                obj.A = length(agents);
            end
            
            obj.expSettings = expSettings;
            obj.clusterMembership = clusterMembership;
            obj.compareIdx = compareIdx;
            
            obj.scrsz = get(0,'ScreenSize');            
            
            import graphs.*
            
            if exist('clusterMembership','var') && ~isempty(obj.clusterMembership)
                obj.clustered = true;
            else
                obj.clustered = false;
            end       
            
            if exist('compareIdx', 'var') && compareIdx > 0
                obj.compare = true;
            else
                obj.compare = false;
            end     
                        
            obj.N = obj.expSettings.nSamples;
            obj.Si = obj.expSettings.nInfSensors;
            obj.K = obj.expSettings.K;
            obj.S = obj.expSettings.nSensors();       
            
            obj.clusterChanges = {obj.A};            
            
            if ~exist('legendStrings','var') || isempty(legendStrings)
                
                obj.agentStrings = {obj.A};
                for a=1:obj.A
                    obj.agentStrings{a} = sprintf('agent %i', a);        
                end                
                
                if obj.clustered
                    obj.clusterStrings = {obj.K*2};
                    for k=1:obj.K
                        obj.clusterStrings{k*2-1} = sprintf('agent in cluster %i', k);        
                        obj.clusterStrings{k*2} = sprintf('change to cluster %i', k);
                    end                    
                end
            else
                obj.agentStrings = legendStrings;
            end
            

            if obj.clustered
                obj.colourSet = graphs.createColourSet(expSettings.K);
            else
                obj.colourSet = graphs.createColourSet(length(obj.agentStrings));
            end             
            
            
            if exist('changes', 'var') && ~isempty(changes)
                obj.changes = changes;
            else
                obj.changes = {[] [] [] []};
            end
        end
                
        function drawGraphs(obj, basePost, baseLogodds, correctLabels, windowLabel)

            obj.drawPostGraph(basePost, windowLabel);
            obj.drawLogoddsGraph(baseLogodds, windowLabel);
            obj.drawErrorGraph(correctLabels, basePost, windowLabel);
            
            if isempty(obj.agents) == false
                obj.drawSensorSubscriptions();
            end
            
            if obj.clustered
               obj.drawClusterChanges(); 
            end
        end
        
        function selectAndPlot(obj, agentIdx, results, comparison)
            if obj.clustered
                obj.clusterChanges{agentIdx} = graphs.plotResultsClusterColours(...
                    agentIdx, obj.clusterMembership, results, obj.colourSet, obj.N, comparison);            
            else
                graphs.plotBasic(agentIdx, results, obj.colourSet, size(results, 2), 2, comparison);
            end
        end     
 
        function drawResultHistogram(obj, resultsPost, windowLabel, dataLegend)
            axisDim = [0 size(resultsPost, 2) ...
                min(reshape(resultsPost, numel(resultsPost), 1))-0.2...
                max(reshape(resultsPost, numel(resultsPost), 1))+0.2];
            
            nBins = 4; 
            centres = [-0.5 0.5 1.5 2.5];
            counts = hist(resultsPost(1,:), centres);
                        
            for r=3
                nSubBins = max(resultsPost(r,:))-min(resultsPost(r,:)) + 1;
                
                if nSubBins > 10
                    nSubBins = 10;
                end
                interval = (max(resultsPost(r,:)) - min(resultsPost(r,:))) / nSubBins;
                subCentres = (double(1:int32(nSubBins)) .* interval) + min(resultsPost(r,:));
                subCentres = subCentres - (interval/2);
                
                if strcmp(obj.agentStrings{r}, 'Decision Count')
                    subCentres = [10 20 30 50];
                end
                
                subCounts = zeros(nBins, length(subCentres));                
                
                prevStart = 1;
                for c=1:nBins
                    subCounts(c, :) = hist(resultsPost(r, prevStart:counts(c)), subCentres);
                end                
                
                figure;
                bar(centres, subCounts, 1, 'stack');
                               
                legend(dataLegend);
                xlabel([obj.agentStrings(1) ' score']);
                ylabel('No. Assets');                
            end
        end        
        
        function drawPostGraph(obj, resultsPost, windowLabel)
            axisDim = [0 size(resultsPost, 2) ...
                min(reshape(resultsPost, numel(resultsPost), 1))-0.2...
                max(reshape(resultsPost, numel(resultsPost), 1))+0.2];
            obj.drawPerfGraph(resultsPost, 'p(c1|x)', axisDim, 'post', windowLabel);
        end

        function drawLogoddsGraph(obj, resultsLogodds, windowLabel)
            if strcmp(obj.expSettings.agentType, 'dlr')
                logit_min = -10;
                logit_max = 10;
            else
                logit_min = -500;
                logit_max = 500;
            end
            axisDim = [0 size(resultsLogodds, 2) logit_min logit_max];
            
            obj.drawPerfGraph(resultsLogodds, 'a', axisDim, 'logodds', windowLabel);
        end
        
        function drawErrorGraph(obj, correctLabels, resultsPost, windowLabel)
            %ERROR OF EACH AGENT, SHOWS CLUSTER MEMBERSHIP, INFORMEDNESS
            
            %fix unknown classification labels - just make them 0.5 for
            %now.
            correctLabels(correctLabels==0, :) = 0.5;
            
            %we don't use the value stored in exp settings for num agents,
            %as the number of "agents" here may refer to the number of 
            %combinations of real agents.
            obj.A = size(resultsPost, 1);         
            errors = obj.calculateBayesErrors(resultsPost, size(resultsPost, 2), correctLabels);
            
            axisDim = [0 size(resultsPost, 2) -0.1 1.1];
            
            obj.drawPerfGraph(errors, 'Absolute Error', axisDim, 'absError', windowLabel);
        end
        
        function drawPerfGraph(obj, results, titleLabel, axisDim, graphTypeLabel, windowLabel)
            %POSTERIOR CLASS PROBS FOR EACH AGENT, SHOWS CLUSTER MEMBERSHIP, INFORMEDNESS
            obj.A = size(results, 1);         
            
            fs = 16; %font size
            
            obj.A = size(results, 1);
            
            %h = figure('Position',[1 obj.scrsz(4) obj.scrsz(3) obj.scrsz(4)],...
            h=figure('Position', [1 1 obj.scrsz(3) obj.scrsz(4)-100],...
                'Visible', obj.expSettings.graphVisible);
            
            for a=1:obj.A 

                %skip the data that is just for comparison, not a real result, if
                %it is set
                if obj.compare && a == obj.compareIdx
                    continue
                end

                label = sprintf('%s for %s', titleLabel, obj.agentStrings{a});

                %select a suitable subplot layout
                if obj.expSettings.saveImgs
                    if obj.compare
                        graphIdx = a;
                        if a > obj.compareIdx
                            graphIdx = a-1;
                        end                
                        subplot(obj.A-1, 1, graphIdx);
                    else
                        subplot(obj.A, 1, a);
                    end
                elseif obj.compare
                    graphIdx = a;
                    if a > obj.compareIdx
                        graphIdx = a-1;
                    end
                    subplot( floor((obj.A)/2), 2, graphIdx);
                else
                    subplot( floor((obj.A+1)/2), 2, a);
                end

                %select and plot a graph
                if obj.compare
                    obj.selectAndPlot(obj.compareIdx, results, true);
                end
                set(gca, 'FontSize', fs);
                obj.selectAndPlot(a, results, false);
                set(gca, 'FontSize', fs);
                if ~isempty(obj.agents)
                    agent = obj.agents(1, a);
                    agent = agent{1};
                    graphs.plotInformedness(agent, obj.Si, obj.N, ...
                        obj.changes, 0, 1);
                else
                    %if we can't plot the informedness of the agents we can
                    %still show the sensor changepoints.
                    graphs.plotSensorChangepoints(obj.changes, true);
                    set(gca, 'FontSize', fs);
                end

                if exist('obj.expSettings.corruptLabels', 'var');
                    agent_corruptLabels = obj.expSettings.corruptLabels(find(obj.expSettings.corruptLabels(:, 2)==a), 1:2);
                    plot(agent_corruptLabels(:, 1), ones(length(agent_corruptLabels),1), 'x', ...
                        'MarkerSize', 5, 'Color', [1 0 1]);

                    agent_missingLabels = obj.expSettings.missingLabels(find(obj.expSettings.missingLabels(:, 3)==a), 1:2);
                    plot(agent_missingLabels(:, 1), agent_missingLabels(:, 2), '*', ...
                        'MarkerSize', 5, 'Color', [0 1 1]);

                    agent_flippedLabels = find(obj.expSettings.flippedLabels(:, a));
                    plot(agent_flippedLabels(:, 1), ones(length(agent_flippedLabels)), 'o', ...
                        'Markersize', 6, 'Color', [0.8 0.8 0.2]);
                end
                hold on
                if obj.N <= 250
                    set(gca, 'XTick', 0:10:obj.N);
                elseif obj.N <= 500
                    set(gca, 'XTick', 0:20:obj.N);
                end
                axis(axisDim);
                set(gca, 'FontSize', fs);
                title(label);    
                set(gca, 'FontSize', fs);
                
                if obj.clustered
                    legend(char(obj.clusterStrings));
                else
                    if obj.compare
                        legend(char({obj.agentStrings{obj.compareIdx} obj.agentStrings{a} }));%'changepoints'}));
                    else
                        legend(char({obj.agentStrings{a}}));
                    end
                end
                set(gca, 'FontSize', fs);
                
                hold all
            end

            if obj.expSettings.saveImgs
                
                if ~exist('windowLabel', 'var')
                    windowLabel = 'general';
                end
                
                saveas(h, sprintf('%s/%s-%s-%s', obj.expSettings.topSaveDir, ...
                    obj.expSettings.expLabel, ...
                    windowLabel, ...
                    graphTypeLabel), 'png');
            end
            
             
        end
        
   
        
        function drawSensorSubscriptions(obj)

            % SENSOR SUBSCRIPTION TABLE
            
            hold off
            if ~isempty(obj.agents)
                obj.A = length(obj.agents);%obj.expSettings.nAgents;
                
                h = figure('Visible', obj.expSettings.graphVisible);
                title('Sensor subscriptions by agent');

                for a=1:obj.A

                    X = zeros(1, obj.agents{a}.num_active_sensors);
                    X = X + a;
                    Y = zeros(1, obj.agents{a}.num_active_sensors);

                    sub = obj.agents{a}.s;%subscriptions{a};
                    sensors_counted = 0;
                    for s=1:obj.S
                        if sub(1, s) == 1
                            sensors_counted = sensors_counted + 1;
                            Y(1, sensors_counted) = s;
                        end
                    end

                    plot(X, Y,...
                    'o', ...
                    'LineWidth',10, ...
                    'MarkerSize',10);

                    hold all
                end
                legend(char(obj.agentStrings));
                if obj.expSettings.saveImgs
                    saveas(h, sprintf('%s/%s-subscriptions', obj.expSettings.topSaveDir, obj.expSettings.expLabel), 'png');
                end
            end
        end
        
        function drawClusterChanges(obj)

            %AGENT CLUSTER CHANGES

            obj.A = length(obj.agents);

            hold off 
            %h = figure('Position',[1 obj.scrsz(4) obj.scrsz(3) obj.scrsz(4)], ...
            h=figure('Position', [1 1 obj.scrsz(3) obj.scrsz(4)],...
                'Visible', obj.expSettings.graphVisible);

            subplot(2, 1, 1);

            cluster_changepoints_y = zeros(obj.N, 1);

            obj.colourSet = graphs.createColourSet(obj.A);

            y_max = 0;

            for a=1:obj.A

                X = obj.clusterChanges{a};

                changes_y_a = zeros(obj.N, 1);

                for i=1:length(X)
                    cluster_changepoints_y(X{i}, 1) = cluster_changepoints_y(X{i}, 1) + 1;

                    if cluster_changepoints_y(X{i}, 1) > y_max
                        y_max = cluster_changepoints_y(X{i}, 1);
                    end

                    changes_y_a(X{i}, 1) = 1;
                end
                
                obj.colourSet(a, :)
                
                plot(1:obj.N, changes_y_a, 'Color', obj.colourSet(a, :));
                hold all;

            end
            axis([0 obj.N -0.5 1.5]);
            
            legend(char(obj.agentStrings));
            title('Cluster changes by agent');

            subplot(2, 1, 2);
            plot(1:obj.N, cluster_changepoints_y);
            axis([0 obj.N -0.5 obj.S]);
            title('Sum of cluster changes');

            hold all
        
            if ~isempty(obj.changes)
                graphs.plotSensorChangepoints(obj.changes);
                if obj.expSettings.saveImgs
                    saveas(h, sprintf('%s/%s-changes', ...
                        obj.expSettings.topSaveDir, ...
                        obj.expSettings.expLabel), 'png');
                end
            end    
        end

            %{
            %AGENT PREDICTION ERROR ALL ON ONE GRAPH - bit hard to read!

            %showing graphs for last test data points only.
            hold off
            h = figure;

            obj.colourSet = graphs.createColourSet(obj.A);

            for a=1:obj.A
                plot(1:obj.N, abs(transpose(correct_labels(1:obj.N, 1)) - basePost(a, 1:obj.N)),...
                    '.-', ...
                    'LineWidth', 1, ...
                    'MarkerSize', 20, ...
                    'Color', obj.colourSet(a, :) ...
                );

                title('Prediction error for all agents');    
                hold all
            end

            %}
    end
end

