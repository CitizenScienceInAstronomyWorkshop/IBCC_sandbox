classdef ClusterGraph < handle
    
    properties
        expSettings
        
        clusterer
        
        scrsz
        
        A
        N
        K
        
        sensorChanges
        
        colourSet
    end
    
    methods
        
        function obj = ClusterGraph(expSettings, changes, clusterer)
            obj.expSettings = expSettings;
            obj.clusterer = clusterer;
            obj.scrsz = get(0,'ScreenSize');
            
            obj.A = obj.expSettings.nAgents;
            obj.N = obj.expSettings.nSamples;
            obj.K = obj.expSettings.K;

            obj.sensorChanges = changes;
            
            obj.colourSet = createColourSet(obj.K);
            
            import graphs.*;
        end

        function drawGraphs(obj, resultsPost)
            obj.drawClusterMembership();
            obj.drawAgentPerfGraph(resultsPost);
            obj.drawClusterStats();
        end
        
        function h = createHalfScreenFigure(obj)
            %h = figure('Position',[1 obj.scrsz(4) obj.scrsz(3)/2 obj.scrsz(4)],...
            %            'Visible', obj.expSettings.graphVisible); 
            h = figure('Position',[1 1 obj.scrsz(3)/2 obj.scrsz(4)],...
                        'Visible', obj.expSettings.graphVisible);             
        end
    
        function drawClusterMembership(obj)
            %CLUSTER MEMBERSHIP OF ALL AGENTS

            %plot the agent membership. Spread the cluster membership values out so we
            %can see the separate lines
            h = obj.createHalfScreenFigure();
            subplot(2, 1, 1);

            legendStrings ={obj.A};

            for a=1:obj.A
                plot(1:obj.N, obj.clusterer.agentClusterMembership(a, :)+(a/(obj.A*4)), '.-', 'MarkerSize',5, 'LineWidth', 1);
                hold all
                legendStrings{a} = sprintf('agent %i', a);
            end
            legend(char(legendStrings));
            title('cluster membership of agents over n steps');

            subplot(2, 1, 2);
            graphs.plotSensorChangepoints(obj.sensorChanges);

            if exist('corruptLabels','var')
                plot(obj.expSettings.corruptLabels(:, 1), obj.expSettings.corruptLabels(:, 2), 'x', ...
                    'MarkerSize', 5, 'Color', [1 0 1]);
                hold all
                plot(obj.expSettings.missingLabels(:, 1), obj.expSettings.missingLabels(:, 2), '*', ...
                    'MarkerSize', 5, 'Color', [0 1 1]);
            end
            title('changepoints marked shown as black upward/downward triangles to indicate sensor  being switched on/off respectively');


            if obj.expSettings.saveImgs
                saveas(h, sprintf('%s/%s-membership', obj.expSettings.topSaveDir,...
                    obj.expSettings.expLabel), 'png');
            end            
        end
        
        
        
        function drawAgentPerfGraph(obj, resultsPost)
            %RESULTS OF ALL AGENTS, COLOURED BY CLUSTER MEMBERSHIP, SHOWING
            %CHANGEPOINTS
            h = obj.createHalfScreenFigure();
            subplot(2, 1, 1);

            for a=1:obj.A
                graphs.plotResultsClusterColours(a, obj.clusterer.agentClusterMembership, ...
                    resultsPost, obj.colourSet, obj.N);
                axis([0 obj.N -1 2 ]);
            end

            title('p(c1|x) of all agents, coloured by cluster');

            subplot(2, 1, 2);
            %draw the changepoint on the graph as well
            graphs.plotSensorChangepoints(obj.sensorChanges);

            if exist('corruptLabels','var')
                plot(obj.expSettings.corruptLabels(:, 1), obj.expSettings.corruptLabels(:, 2), 'x', ...
                    'MarkerSize', 5, 'Color', [1 0 1]);
                plot(obj.expSettings.missingLabels(:, 1), obj.expSettings.missingLabels(:, 2), '*', ...
                    'MarkerSize', 5, 'Color', [0 1 1]);
            end
            title('changepoints marked shown as black upward/downward triangles to indicate sensor  being switched on/off respectively');

            if obj.expSettings.saveImgs
                saveas(h, sprintf('%s/%s-postall', obj.expSettings.topSaveDir, obj.expSettings.expLabel), 'png');
            end
            
        end
        
        
        function drawClusterStats(obj)
            %CLUSTER SIZES, ERRORS, CONFIDENCE, MEANS
            h = obj.createHalfScreenFigure();

            legendStrings = {obj.K};
            for k=1:obj.K
                legendStrings{k} = sprintf('cluster %i', k);
            end

            %plot the error of each cluster over n iterations
            for k=1:obj.K
                colour = obj.colourSet(k, :);

                subplot(4, 1, 1);
                plot(1:obj.N, obj.clusterer.cluster_errors(k, 1:obj.N), ...
                    'LineWidth', 1, 'Color', colour);

                title('error of clusters');   

                hold all

                subplot(4, 1, 3);
                plot(1:obj.N, obj.clusterer.cluster_sizes(k, 1:obj.N), ...
                    'LineWidth', 1, 'Color', colour);

                title('sizes of clusters');

                hold all

                subplot(4, 1, 2);
                plot(1:obj.N, obj.clusterer.cluster_confidences(k, 1:obj.N), ...
                    'LineWidth', 1, 'Color', colour);

                title('confidence in clusters');

                hold all

                subplot(4, 1, 4);
                plot(1:obj.N, obj.clusterer.cluster_means(k, 1:obj.N), 'o', ...
                    'LineWidth', 1, 'Color', colour, 'MarkerFaceColor', colour);

                title(sprintf('means of each cluster (mean %s values)', obj.expSettings.clusterData));

                if strcmp(obj.expSettings.clusterData, 'logodds')
                    axis([0 obj.N -500 500]);
                elseif strcmp(obj.expSettings.clusterData, 'post')
                    axis([0 obj.N -0.1 1.1]);
                end

                hold all

                legend(char(legendStrings));
            end

            if obj.expSettings.saveImgs
                saveas(h, sprintf('%s/%s-clusterstats', obj.expSettings.topSaveDir, obj.expSettings.expLabel), 'png');
            end

            hold off            
            
        end



        %{
        h = figure;

        legendStrings = {};
        for s=1:nSensorsC
            legendStrings{s} = sprintf('sensor %i', s);
        end

        colourSet = createColourSet(nSensors);
        for k=1:K        
            sensor_acc_credit = scoring.calculate_cluster_sensor_credit(agents, A, nSensors, N, classifications, basePost, clusters, k);
            subplot(K, 1, k);

            for s=1:nSensors

                plot( 1:N, sensor_acc_credit(s, :), 'Color', colourSet(s,:) );
                title(sprintf('Accumulating Sensor Credit for cluster %d', k));
                hold all
            end
        end

        legend( char(legendStrings) );

        hold off
        %}
    end
end