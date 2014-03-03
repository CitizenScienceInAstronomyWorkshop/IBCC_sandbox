classdef VariableSettingGraph
    %VARIABLESETTINGGRAPH Graph showing level of error when a setting for a
    %variable is altered.
    %   Can take as input the results of a number of runs with the same
    %   setting and with different settings, and plots the error level as
    %   the variable changes.
    
    properties
        label = 'error rates with different variable settings';
        nSamples
        colourSet
        
        fontsize = 24;
        errorWindowSize = 1;
        
        labelXAxis = true;
        labelYAxis = true;
        
        yLimit = [];
        xLimit = [];
        
        shade = false;
        dashed = false;
    end
    
    methods
        function obj = VariableSettingGraph(label, nSamples, nCombiners)
            obj.label = label;
            obj.nSamples = nSamples;
            if nargin > 2
                nColours = nCombiners;
                if nColours < 9
                    nColours = 9;
                end
                obj.colourSet = graphs.createColourSet(nColours);
                    
                obj.colourSet(4, :) = [1 0.7 0];
                %obj.colourSet(5, :) = [0.9 0.9 0];
                obj.colourSet(6, :) = [0.1 0.1 1];
                
                %obj.colourSet(8, :) = [0 0.8 0];
            end
        end
        
        function drawParzenGraph(obj, useMeans, combinedPosts, allLabels, nRepeats, ...
                varSettings, parzenWindowSize, varSettingsName, combIdx, ...
                reuseGraph, legendStrings)        
            
            if nargin < 6
                parzenWindowSize = 0.002;
            end
            
            if nargin < 8
                combIdx = 1;
            end
            
            nSettings = size(combinedPosts, 2);
            nDatasets = size(combinedPosts, 1);
            nDataSamples = size(allLabels, 2);
            
            if useMeans
                nPoints = nDatasets*nRepeats;
            else
                nPoints = nDatasets*nRepeats*nDataSamples;
                
                startIdx = 1;
                lastIdx = nRepeats*nDataSamples;
            end
            
            totalErrors = zeros(nSettings, nPoints);
                        
            for s=1:nSettings
                
                %construct a list of points for each setting
                for d=1:nDatasets
                                        
                    repeats = combinedPosts{d, s};
                    repeats = repeats(:, :, combIdx);
                    
                    display('assume variable settings give number of known targets');
                    
                    nKnowns = varSettings(s);
                    knownLabels = ExpRunner.getKnownTargets(allLabels(d, :), nKnowns, false);
                    
                    labels = allLabels(d, 1:size(repeats, 2));
                    labels(knownLabels~=-1) = -1;
                    
                    labels = ones(nRepeats, 1)*labels;
                    errors = abs(repeats - labels);
                    errors(find(labels==-1)) = 0;
                    
                    totalSamples = sum(labels~=-1, 2);
                    
                    %take the mean
                    if useMeans
                        errors = sum(errors, 2) ./ totalSamples;             
                                        
                        totalErrors(s, (d-1)*nRepeats+1:d*nRepeats) = ...
                            errors';
                    else
                        totalErrors(s, startIdx:lastIdx) = reshape(errors, 1, nDataSamples*nRepeats);
                        startIdx = lastIdx + 1;
                        lastIdx = lastIdx + nDataSamples*nRepeats;
                    end
                end
            end        
            
            if nargin < 11
                legendStrings = [];
            end
            
            if nargin > 9
                obj.drawParzenGivenError(totalErrors, [], parzenWindowSize, reuseGraph, ...
                    varSettings, varSettingsName, combIdx, legendStrings);
            else
                obj.drawParzenGivenError(totalErrors, [], parzenWindowSize, false, ...
                    varSettings, varSettingsName, combIdx, legendStrings);
            end
        end
        
        function drawParzenGivenError(obj, totalErrors, errorVars, parzenWindowSize, ...
                reuseGraph, varSettings, varSettingsName, combIdx, legendStrings)             
            
            %[M, C, P] = graphs.parzen(totalErrors');
            
            nGraphPoints = 101;
            nSettings = size(totalErrors, 1);
            nPoints = size(totalErrors, 2);
            
            M = totalErrors;
            if isempty(errorVars)
                C = ones(nSettings, nPoints) .* parzenWindowSize;
            else
                C = errorVars .*parzenWindowSize;
            end
            
            if nargin >= 4 && reuseGraph
                hold all;
            else
                h = figure;
                hold off;
            end
            
            %calculate function f(x), which is graph height dependent on
            %density of samples with error value close to x.

            X = obj.errorWindowSize * (0:(nGraphPoints-1)) ./ (nGraphPoints-1);
            plotData = zeros(nGraphPoints, nSettings);
            
            for s=1:nSettings
                
                contributions = normpdf(ones(nPoints, 1)*X,...    
                    M(s, :)' * ones(1,nGraphPoints),...
                    C(s, :)' * ones(1, nGraphPoints));
                
                
                fX = sum(contributions, 1) ./ nPoints;
                
                plotData(1:nGraphPoints, s) = fX;
                %plot(X, fX);
                %plot3(varSettings(s).*ones(1, nGraphPoints), X, fX);
                %hold on
            end
            
            fontsize = obj.fontsize;
            
            if length(varSettings) == 1
                if combIdx <= length(obj.colourSet)
                    col = obj.colourSet(combIdx, :);
                end
                linestyle = '-';
                linewidth = 3;
                marker = 'none';
%                 switch combIdx 
%                     case 1
%                         marker = 'x';
%                     case 2
%                         marker = '*';
%                     case 3
%                         marker = '^';
%                     case 5 
%                         marker = 's';
%                         linestyle = '--';
%                     case 6
%                         marker = 'o';
%                         linestyle = '--';
%                     case 8
%                         marker = '+';
%                         linestyle = '--';
%                     otherwise
%                         marker = '';
%                 end
                                
                if obj.shade == false
                    plot(X, plotData, 'Color', col, 'LineWidth', linewidth, ...
                        'LineStyle', linestyle, 'Marker', marker, 'MarkerSize', 6);
                    
                else
                    area(X, plotData);
                end
                             
                set(gca, 'FontSize', fontsize);
                
                if nargin > 8 && ~isempty(legendStrings)
                    legend(legendStrings);
                end
                if nargin < 9 || ~isempty(legendStrings) || ~reuseGraph
                    if obj.labelXAxis
                        xlabel('Mean absolute error for all data points', 'fontsize', fontsize);
                    else
                        %set(gca,'XTickLabel',{''})
                    end
                end

                if obj.labelYAxis
                    ylabel('Parzen density', 'fontsize', fontsize); 
                end


                
                set(gca,'ZTickLabel',{''})
                
                if ~isempty(obj.yLimit)
                    set(gca, 'YLim', obj.yLimit);
                end
                if ~isempty(obj.xLimit)
                    set(gca, 'XLim', obj.xLimit);
                end                
                box on
            else
                
                surf(varSettings, X, plotData, 'FaceColor', 'interp', 'EdgeColor', 'interp', ...
                    'EdgeAlpha', 0.7);%'MeshStyle', 'both', 

                set(gca, 'FontSize', fontsize);
                set(gca, 'YTick', [0.2 0.4 0.5 0.6 0.8 1.0]);
                %set(gca, 'XMinorTick', 'on');
                %set(gca, 'YMinorTick', 'on');
                %set(gca,'ZTickLabel',{''})
                if ~isempty(obj.xLimit)
                    set(gca, 'XLim', obj.xLimit);
                end  
                if ~isempty(obj.yLimit)
                    set(gca, 'YLim', obj.yLimit);
                end                
                hold on;
                
                if nargin < 9 || ~reuseGraph
                    ylabel('Mean absolute error', 'fontsize', fontsize);
                    xlabel(varSettingsName, 'fontsize', fontsize);
                    zlabel('Parzen density', 'fontsize', fontsize);
                end
            
                view([-120 50]);
            end
            if ~isempty(obj.label)
                title(obj.label, 'fontsize', fontsize);
            end
        end        
        
        function drawErrorHist(obj, totalErrors, intervalSize, reuseGraph)             
                        
            nBins = obj.errorWindowSize / intervalSize;
            
            nSettings = size(totalErrors, 1);
            nPoints = size(totalErrors, 2);
            
            M = totalErrors;
            
            if nargin >= 4 && reuseGraph
                hold all;
            else
                figure;
                hold off;
            end

            hist(totalErrors, nBins);

            set(gca, 'FontSize', obj.fontsize);

            if nargin > 8 && ~isempty(legendStrings)
                legend(legendStrings);
            end
            if nargin < 9 || ~isempty(legendStrings) || ~reuseGraph
                if obj.labelXAxis
                    xlabel('Mean Absolute Error (Bayes Error)', 'fontsize', obj.fontsize);
                else
                    %set(gca,'XTickLabel',{''})
                end
            end

            if obj.labelYAxis
                ylabel('Number of Base Classifiers', 'fontsize', obj.fontsize); 
            end

            set(gca,'ZTickLabel',{''})

            if ~isempty(obj.yLimit)
                set(gca, 'YLim', obj.yLimit);
            end
            if ~isempty(obj.xLimit)
                set(gca, 'XLim', obj.xLimit);
            end                
            box on
            if ~isempty(obj.label)
                title(obj.label, 'fontsize', obj.fontsize);
            end
        end             
        
        function [combDatasetsMean, combDatasetsDev, combDatasetsMin, combDatasetsMax] ...
                = error2dRuns(obj, combinedPosts, allLabels, acrossDatasets) 
            
            nSettings = size(combinedPosts, 2);
            nDatasets = size(combinedPosts, 1);
            varRepeats = zeros(nDatasets, nSettings);
            meanRepeats = zeros(nDatasets, nSettings);
            minRepeats = zeros(nDatasets, nSettings);
            maxRepeats = zeros(nDatasets, nSettings);
            
            totalRuns = zeros(1, nSettings);
            for s=1:nSettings
                
                for d=1:nDatasets
                    
                    
%                     if d==2
%                         display('skipping dataset 2 - no informed agents');
%                         continue
%                     end
                    
                    if isempty(combinedPosts{d, s})    
                        continue
                    end
                    
                    repeats = combinedPosts{d, s};
                    nRepeats = size(repeats, 1);
                    
                    labels = ones(nRepeats, 1)*allLabels(d, :);
                    errors = abs(repeats - labels);
                                        
                    totalRuns(1, s) = totalRuns(1, s) + nRepeats;
                    
                    meanErrors = sum(errors, 1) ./ nRepeats;
                    varErrors = sum((errors - ones(nRepeats, 1)*meanErrors) .^ 2, 1);
                    varErrors = sum(varErrors, 2) ./ (nRepeats*obj.nSamples);
                    meanErrors = sum(meanErrors, 2) ./ obj.nSamples;
                    
                    totalErrors = sum(errors, 2) ./ obj.nSamples;
                    minRepeats(d, s) = min(totalErrors);
                    maxRepeats(d, s) = max(totalErrors);
                         
                    %sum across each target label
                    if acrossDatasets
                        varRepeats(d, s) = varErrors * nRepeats;
                        meanRepeats(d, s) = meanErrors * nRepeats;
                    else
                        varRepeats(d, s) = varErrors;
                        meanRepeats(d, s) = meanErrors;
                    end
                end
            end               
            
            %combine across datasets
            if acrossDatasets
                combDatasetsDev = sum(varRepeats, 1) ./ totalRuns;
                combDatasetsMean = sum(meanRepeats, 1) ./ totalRuns;
                
                combDatasetsMin = min(minRepeats);
                combDatasetsMax = max(maxRepeats);
            else
                combDatasetsDev = varRepeats;
                combDatasetsMean = meanRepeats;
                
                combDatasetsMin = minRepeats;
                combDatasetsMax = maxRepeats;
            end
            
            combDatasetsDev = combDatasetsDev .^ 0.5;
        end 
        
        function [comb1Mean, comb1Dev, comb1Min, comb1Max, ...
                comb2Mean, comb2Dev, comb2Min, comb2Max] ...
                = error3dRuns(obj, combinedPosts, allLabels, acrossDatasets) 
            
            nSettings1 = size(combinedPosts, 2);
            nSettings2 = size(combinedPosts, 3);
            nDatasets = size(combinedPosts, 1);
            varRepeats = zeros(nSettings1, nSettings2, nDatasets);
            meanRepeats = zeros(nSettings1, nSettings2, nDatasets);
            minRepeats = zeros(nSettings1, nSettings2, nDatasets);
            maxRepeats = zeros(nSettings1, nSettings2, nDatasets);
            
            totalRuns = zeros(nSettings1, nSettings2);
            for s=1:nSettings1
                for t=1:nSettings2
                    for d=1:nDatasets
                        
                        if d==1
                            display('skipping dataset 2 for lambda comparisons - no informed agents');
                            continue;
                        end
                        
                        if isempty(combinedPosts{d, s})    
                            continue
                        end

                        repeats = combinedPosts{d, s, t};
                        nRepeats = size(repeats, 1);
                        labels = ones(nRepeats, 1)*allLabels(d, :);
                        errors = abs(repeats - labels);

                        totalRuns(s, t) = totalRuns(s, t) + nRepeats;

                        meanErrors = sum(errors, 1) ./ nRepeats;
                        varErrors = sum((errors - ones(nRepeats, 1)*meanErrors) .^ 2, 1);
                        varErrors = sum(varErrors, 2) ./ (nRepeats*obj.nSamples);
                        meanErrors = sum(meanErrors, 2) ./ obj.nSamples;

                        totalErrors = sum(errors, 2) ./ obj.nSamples;
                        
                        minRepeats(s, t, d) = min(totalErrors);
                        maxRepeats(s, t, d) = max(totalErrors);
                    
                        %sum across each target label
                        if acrossDatasets
                            varRepeats(s, t, d) = varErrors * nRepeats;
                            meanRepeats(s, t, d) = meanErrors * nRepeats;
                        else
                            varRepeats(s, t, d) = varErrors;
                            meanRepeats(s, t, d) = meanErrors;
                        end
                    end
                end
            end               

            %combine across datasets
            if acrossDatasets
                comb1Dev = zeros(nSettings2, nSettings1);
                comb1Mean = zeros(nSettings2, nSettings1);
                comb1Min = zeros(nSettings2, nSettings1);
                comb1Max = zeros(nSettings2, nSettings1);
                
                for t=1:nSettings2
                    comb1Dev(t, :) = sum(varRepeats(:, t, :), 3)' ./ totalRuns(:, t)';
                    comb1Mean(t, :) = sum(meanRepeats(:, t, :), 3)' ./ totalRuns(:, t)';
                    comb1Min(t, :) = min(minRepeats(:, t, :), [], 3);
                    comb1Max(t, :) = max(maxRepeats(:, t, :), [], 3);
                end
                
                comb2Dev = zeros(nSettings1, nSettings2);                
                comb2Mean = zeros(nSettings1, nSettings2);
                comb2Min = zeros(nSettings1, nSettings2);
                comb2Max = zeros(nSettings1, nSettings2);
                
                for s=1:nSettings1
                    comb2Dev(s, :) = sum(varRepeats(s, :, :), 3) ./ totalRuns(s, :);
                    comb2Mean(s, :) = sum(meanRepeats(s, :, :), 3) ./ totalRuns(s, :);
                    comb2Min(s, :) = min(minRepeats(s, :, :), [], 3);
                    comb2Max(s, :) = max(maxRepeats(s, :, :), [], 3);
                end                
            else
                comb1Dev = zeros(nSettings2*nDatasets, nSettings1);
                comb1Mean = zeros(nSettings2*nDatasets, nSettings1);
                comb1Min = zeros(nSettings2*nDatasets, nSettings1);
                comb1Max = zeros(nSettings2*nDatasets, nSettings1);
                
                for d=1:nDatasets
                    rows = ((d-1) .* (0:nSettings2-1)) + 1;
                    comb1Dev(rows, :) = varRepeats(:, :, d)';
                    comb1Mean(rows, :) = meanRepeats(:, :, d)';
                    comb1Min(rows, :) = minRepeats(:, :, d)';
                    comb1Max(rows, :) = maxRepeats(:, :, d)';                    
                end
                
                comb2Dev = zeros(nSettings1*nDatasets, nSettings2);
                comb2Mean = zeros(nSettings1*nDatasets, nSettings2);
                comb2Min = zeros(nSettings1*nDatasets, nSettings2);
                comb2Max = zeros(nSettings1*nDatasets, nSettings2);
                for d=1:nDatasets
                    rows = ((d-1) .* (0:nSettings2-1)) + 1;
                    comb2Dev(rows, :) = varRepeats(:, :, d);
                    comb2Mean(rows, :) = meanRepeats(:, :, d);
                    comb2Min(rows, :) = minRepeats(:, :, d);
                    comb2Max(rows, :) = maxRepeats(:, :, d);
                end
            end
            
            comb1Dev = comb1Dev .^ 0.5;
            comb2Dev = comb2Dev .^ 0.5;            
         end  
        
        function [combDatasets] = disagreement2dRuns(obj, combinedPosts, acrossDatasets) 
            
            nSettings = size(combinedPosts, 2);
            nDatasets = size(combinedPosts, 1);
            varRepeats = zeros(nDatasets, nSettings);
            
            totalRuns = zeros(1, nSettings);
            for s=1:nSettings
                
                for d=1:nDatasets
                    if isempty(combinedPosts{d, s})    
                        continue
                    end
                    
                    repeats = combinedPosts{d, s};
                    nRepeats = size(repeats, 1);
                    
                    totalRuns(1, s) = totalRuns(1, s) + nRepeats;
                    
                    means = sum(repeats, 1) ./ nRepeats;
                    sqDiffs = (repeats - ones(nRepeats, 1)*means) .^ 2;
                    variance = sum(sum(sqDiffs, 1), 2) ./ (nRepeats*obj.nSamples);                    
                                        
                    %sum across each target label
                    if acrossDatasets
                        varRepeats(d, s) = variance * nRepeats;
                    else
                        varRepeats(d, s) = variance;
                    end
                end
            end               
            
            %combine across datasets
            if acrossDatasets
                combDatasets = sum(varRepeats, 1) ./ totalRuns;
            else
                combDatasets = varRepeats;
            end
            
            combDatasets = combDatasets .^ 0.5;
        end   
                
        function [combDatasets1, combDatasets2] = disagreement3dRuns(obj, combinedPosts, acrossDatasets) 
            %combDatasets1 contains a different dataset for each value of
            %setting 2, combDatasets2 contains different result set for
            %each value of setting 1.
            nSettings1 = size(combinedPosts, 2);
            nSettings2 = size(combinedPosts, 3);
            nDatasets = size(combinedPosts, 1);
            varRepeats = zeros(nSettings1, nSettings2, nDatasets);
            
            totalRuns = zeros(nSettings1, nSettings2);
            for s=1:nSettings1
                for t=1:nSettings2
                    for d=1:nDatasets
                        if isempty(combinedPosts{d, s})    
                            continue
                        end

                        repeats = combinedPosts{d, s, t};
                        nRepeats = size(repeats, 1);

                        totalRuns(s, t) = totalRuns(s, t) + nRepeats;

                        means = sum(repeats, 1) ./ nRepeats;
                        sqDiffs = (repeats - ones(nRepeats, 1)*means) .^ 2;
                        variance = sum(sum(sqDiffs, 1), 2) ./ (nRepeats*obj.nSamples);                    

                        %sum across each target label
                        if acrossDatasets
                            varRepeats(s, t, d) = variance * nRepeats;
                        else
                            varRepeats(s, t, d) = variance;
                        end
                    end
                end
            end               
            
            %combine across datasets
            if acrossDatasets
                combDatasets1 = zeros(nSettings2, nSettings1);
                for t=1:nSettings2
                    combDatasets1(t, :) = sum(varRepeats(:, t, :), 3)' ./ totalRuns(:, t)';
                end
                
                combDatasets2 = zeros(nSettings1, nSettings2);
                for s=1:nSettings1
                    combDatasets2(s, :) = sum(varRepeats(s, :, :), 3) ./ totalRuns(s, :);
                end
            else
                combDatasets1 = zeros(nSettings2*nDatasets, nSettings1);
                for d=1:nDatasets
                    rows = ((d-1) .* (0:nSettings2-1)) + 1;
                    combDatasets1(rows, :) = varRepeats(:, :, d)';
                end
                
                combDatasets2 = zeros(nSettings1*nDatasets, nSettings2);
                for d=1:nDatasets
                    rows = ((d-1) .* (0:nSettings2-1)) + 1;
                    combDatasets2(rows, :) = varRepeats(:, :, d);
                end
            end
            
            combDatasets1 = combDatasets1 .^ 0.5;
            combDatasets2 = combDatasets2 .^ 0.5;
        end
        

        function drawDisagreementGraph(obj, combResults, varSettings, settingsForLabels, drawLegend)
            figure;
            plot(varSettings, combResults, 'x-');
            title(obj.label);
            
            if nargin > 3
                legendStrings = {length(settingsForLabels)};
                for i=1:length(settingsForLabels)
                    legendStrings{i} = num2Str(settingsForLabels(i));
                end
                if drawLegend
                    legend(legendStrings);    
                end
            end
        end
        
        function drawErrorGraphMinMax(obj, combResultsMeans, combResultsMin, combResultsMax, ...
                varSettings, settingsForLabels, drawLegend)        
            
            upper = combResultsMax;
            lower = combResultsMin;
            label = sprintf('%s, with min and max', obj.label);
            obj.drawErrorGraph(combResultsMeans, lower, upper, ...
                varSettings, label, settingsForLabels, drawLegend);
            
        end
        
        function drawErrorGraphDev(obj, combResultsMeans, combResultsDev, ...
                varSettings, settingsForLabels, drawLegend)        
            
            upper = combResultsMeans + combResultsDev;
            lower = combResultsMeans - combResultsDev;
            label = sprintf('%s, with standard deviations', obj.label);
            obj.drawErrorGraph(combResultsMeans, lower, upper, ...
                varSettings, label, settingsForLabels, drawLegend);
            
        end
        
        function combinedPosts2d = results3dTo2d(obj, combinedPosts3d, dimToCombine, nRepeats, nSamples)
            %reshape a 3d result grid with results for two sets of settings
            %down to results for one set of settings.
            if dimToCombine == 2
                
%                 combinedPosts2d = reshape(combinedPosts3d, ...
%                     size(combinedPosts3d, 1)*size(combinedPosts3d, 2), ...
%                     size(combinedPosts3d, 3));

                combinedPosts2d = cell(size(combinedPosts3d, 1), size(combinedPosts3d, 3));
                for d=1:size(combinedPosts3d, 1)
                    for s=1:size(combinedPosts3d, 3)
                        
                        repeats = zeros(nRepeats*size(combinedPosts3d, 2), ...
                            nSamples);
                        
                        for t=1:size(combinedPosts3d, 2)
                            repeats( nRepeats*(t-1)+1:nRepeats*t, ...
                                :) = combinedPosts3d{d, t, s};
                        end                        
                        combinedPosts2d{d, s} = repeats;
                    end
                end
                
            elseif dimToCombine == 3
                combinedPosts2d = cell(size(combinedPosts3d, 1), size(combinedPosts3d, 2));
                
                for d=1:size(combinedPosts3d, 1)
                    for s=1:size(combinedPosts3d, 2)
                        
                        repeats = zeros(nRepeats*size(combinedPosts3d, 3), ...
                            nSamples);
                        
                        for t=1:size(combinedPosts3d, 3)
                            repeats( nRepeats*(t-1)+1:nRepeats*t, ...
                                :) = combinedPosts3d{d, s, t};
                        end                        
                        combinedPosts2d{d, s} = repeats;
                    end
                end                
%                 combinedPosts3d = permute(combinedPosts3d, [1, 3, 2]);
%                     
%                 combinedPosts2d = reshape(combinedPosts3d, ...
%                    size(combinedPosts3d, 1)*size(combinedPosts3d, 2), ...
%                    size(combinedPosts3d, 3));       
                %combinedPosts2d = combinedPosts3d(:, :, 1);
            else
                display('cannot combine along this dimension');
            end            
%             
%             if nargin < 4
%                 return;
%             end
%             combLabels = [];
%             
%             for s=1:size(combinedPosts3d, 2)
%                 combLabels = [combLabels; labels];
%             end
        end
                
        function drawErrorGraph(obj, combResultsMeans, lower, upper, ...
                varSettings, label, settingsForLabels, drawLegend)
            colours = graphs.createColourSet(size(combResultsMeans, 1));
            
            %plot the means
            figure;
            for i=1:size(combResultsMeans, 1)
                plot(varSettings, combResultsMeans(i, :), 'x-', 'Color', colours(i, :));
                hold all
            end
            if nargin > 4
                legendStrings = {length(settingsForLabels)};
                for i=1:length(settingsForLabels)
                    legendStrings{i} = num2str(settingsForLabels(i));
                end
                if drawLegend
                    legend(legendStrings);    
                end
            end            
            
            %shade in the deviations or min/max bars
            for i=1:size(combResultsMeans, 1)
                
                filled = [upper(i, :), fliplr(lower(i, :))];
                xPoints = [varSettings, fliplr(varSettings)];

                %plot the data
                fillhandle = fill(xPoints, filled, colours(i, :));
                
                transparency = 0.2;
                %set edge color
                set(fillhandle, 'FaceAlpha', transparency,...
                    'EdgeAlpha', transparency, 'EdgeColor', colours(i, :));
                
                hold all;
                %plot(varSettings, combResultsMeans(i, :)+combResultsDev(i, :), ...
                %    '-', 'Color', colours(i, :));
                %hold all
                %plot(varSettings, combResultsMeans(i, :)-combResultsDev(i, :), ...
                %    '-', 'Color', colours(i, :));
            end
            title(label);
            

        end 
    end
    
    end

