%clear

if ~exist('threeDim', 'var')
    threeDim = true;
end

import settings.*;
variableLength = false;

% vbIbccPaper.exp1; resultsDir = '/homes/49/edwin/matlab/combination/results/vbIbccPaper_temp/exp1';
% vbIbccPaper.exp2; resultsDir = '/homes/49/edwin/matlab/combination/results/vbIbccPaper_temp/exp2';
%vbIbccPaper.exp2inversions; resultsDir = '/homes/49/edwin/matlab/combination/results/vbIbccPaper_temp/exp2inv';
%vbIbccPaper.exp2seq; resultsDir = '/homes/49/edwin/matlab/combination/results/vbIbccPaper/exp2seq';

%bccSettings; expSettings.topSaveDir ='/homes/49/edwin/matlab/combination/data/vbPrsMerged/bcc3'; resultsDir = '/homes/49/edwin/matlab/combination/results/supplement/bcc3';
%bccSettingsHard; expSettings.topSaveDir = '/homes/49/edwin/matlab/combination/data/vbPrsMerged/bcc4Hard'; resultsDir = '/homes/49/edwin/matlab/combination/results/supplement/bcc4Hard';
%bccSettingsExp3; resultsDir = '/homes/49/edwin/matlab/combination/results/supplement/bccExp3'; %expSettings.topSaveDir = '/homes/49/edwin/matlab/combination/data/vbPrsMerged/bccExp3';
% bccSettingsInversions; expSettings.topSaveDir = '/homes/49/edwin/matlab/combination/data/vbPrsMerged/bccInv'; resultsDir = '/homes/49/edwin/matlab/combination/results/supplement/bccInv';
% bccSettingsDynamicsExp3; resultsDir = '/homes/49/edwin/matlab/combination/results/dynamicVB/paper/multi_dynamics'; %variableLength = true; expSettings.topSaveDir = '/homes/49/edwin/matlab/combination/data/dynamicVB/bccDynamicsExp3_2'; 

% bccSettings; resultsDir = '/homes/49/edwin/matlab/combination/results/revisedVB/exp1_test';
%bccSettingsHard; resultsDir = '/homes/49/edwin/matlab/combination/results/revisedVB/exp2';
%bccSettingsExp3; resultsDir = '/homes/49/edwin/matlab/combination/results/revisedVB/exp3'; %expSettings.topSaveDir = '/homes/49/edwin/matlab/combination/data/vbPrsMerged/bccExp3';
%bccSettingsInversions; resultsDir = '/homes/49/edwin/matlab/combination/results/revisedVB/exp4'; 
 
% bccSettingsSimple1; resultsDir = '/homes/49/edwin/matlab/combination/results/simple_base_classifiers/exp1';
%bccSettingsSimple2; resultsDir = '/homes/49/edwin/matlab/combination/results/simple_base_classifiers/exp2';
% bccSettingsSimple3; resultsDir = '/homes/49/edwin/matlab/combination/results/dynamicVB/simple/exp3';
%bccSettingsSimple4; resultsDir = '/homes/49/edwin/matlab/combination/results/simple_base_classifiers/exp4';
% bccSettingsSimpleRandom; resultsDir = '/homes/49/edwin/matlab/combination/results/simple_base_classifiers/expRandom';


if ~exist(resultsDir, 'file')
    mkdir(resultsDir);
end

[labels, combinedPostKnowns, combMethods, basePost] = loadMultTestResults(expSettings);
% combMethods{1} = 'Ibcc-VB';
% combMethods{3} = 'Ibcc-EM';
% combMethods{8} = 'Weighted Sum';
% combMethods{9} = 'Weighted Majority';
% combMethods{10} = 'Ibcc-Gibbs';

labels = labels+1;
noScore = 0;

if length(combMethods) > 9
    nColours = length(combMethods);
else
    nColours = 10;
end
colourSet = graphs.createColourSet(nColours);
scrsz = get(0,'ScreenSize');

% colourSet(2, :) = [0.6 0.45 0.45];
% colourSet(4, :) = [1 0.7 0];
% %obj.colourSet(5, :) = [0.9 0.9 0];
% colourSet(6, :) = [0.2 0 1];
% colourSet(7, :) = [0.4 0.9 0];
% colourSet(8, :) = [0 0.9 1];
% colourSet(9, :) = [0.9 0.3 0.8];
% colourSet(10, :) = [0.5 0 0.7];

%index for weighted majority - use it to flip back stray methods
wmIdx = 3;
% 
% %indexes of proportions of known targets that we want to use
selectedI = [1];

selectedMethods = [1 2 3 4];% 8];%1:length(combMethods);
nMethods = length(selectedMethods);

h = figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)*0.6]);
fs = 16;

selectedMethods(selectedMethods>length(combMethods)) = [];

%draw a separate window for each number of known targets
for i=1:length(selectedI)
    
    p = selectedI(i);
    
        for c=1:length(selectedMethods)
           subplot(length(selectedI), ceil(nMethods), c+(i-1)*nMethods);
            for d=1:expSettings.nDatasets
                datasetOutputs = combinedPostKnowns{d, p};
            
                targets = labels(d, (1:size(datasetOutputs, 2)));                
                
                if exist('dataLengths', 'var')
                    nKnown = expSettings.propKnown(expSettings.iPropKnown) * dataLengths(i);
                else
                    nKnown = expSettings.propKnown(p) * length(targets);
                end
                
                knownLabels = ExpRunner.getKnownTargets(...
                    targets, ...
                    nKnown,...
                    false);
                targets = targets(knownLabels==noScore);
                datasetOutputs = datasetOutputs(:, knownLabels==noScore, :);
                outputs = zeros(expSettings.nRepeats, size(datasetOutputs, 2));
                nRepeats = expSettings.nRepeats;
                
                m = selectedMethods(c);
                if m==10
                    for r=1:nRepeats
                        post = datasetOutputs(r, :, m);
                        row = r;
                                                
                        outputs(row, :) = post;
                        [tpr, fpr, thresholds] = roc(targets-1, outputs(row,:));
                        tpr = [tpr 1];
                        fpr = [fpr 1];    
                        plot(fpr, tpr, 'linewidth', 2, 'clipping', 'off', 'Color', colourSet(mod(d,size(colourSet,1))+1, :));                        
                    end
                    outputs = reshape(outputs', 1, numel(outputs));
                    targets = targets' * ones(1, nRepeats);
                    targets = reshape(targets, 1, numel(outputs));
                    
                else
                    post = datasetOutputs(1, :, m);
                    row = (d-1)*expSettings.nDatasets+1;
                    outputs = post;
                    [tpr, fpr, thresholds] = roc(targets-1, outputs);
                    tpr = [tpr 1];
                    fpr = [fpr 1];                    
                    plot(fpr, tpr, 'linewidth', 2, 'clipping', 'off', 'Color', ...
                        colourSet(mod(d, nColours)+1, :));
                end

                hold all;
            end
            axis([0 1 0 1]);
            if d==1
                legend(combMethods);
                set(gca, 'FontSize', fs);                
            end
            set(gca, 'FontSize', fs);
            xlabel('False Positive Rate');
            set(gca, 'FontSize', fs);
            ylabel('True Positive Rate');
            set(gca, 'FontSize', fs);
            if variableLength == true
%                 nKnown = expSettings.propKnown(p) * expSettings.nSamples;
                title(sprintf('%s with %d known targets', combMethods{m}, nKnown));                        
                set(gca, 'FontSize', fs);        
            else
%                 nKnown = expSettings.propKnown(p) * expSettings.nSamples;
                title(sprintf('%s with %d known targets', combMethods{m}, nKnown));        
                set(gca, 'FontSize', fs);
            end
        end
        
%     figure;
%     for c=1:length(combMethods)
%        subplot(2, ceil(length(combMethods)/2), c);
%        
%        outputs = zeros(expSettings.nRepeats*expSettings.nDatasets, expSettings.nSamples);
%        targets = zeros(expSettings.nRepeats*expSettings.nDatasets, expSettings.nSamples);
%        
%        
%         for d=1:expSettings.nDatasets
%             
%             datasetOutputs = combinedPostKnowns{d, i};
% 
%             nRepeats = expSettings.nRepeats;
%             for r=1:nRepeats
%                 post = datasetOutputs(r, :, c);
%                 row = r + (d-1)*nRepeats;
%                 outputs(row, :) = post;   
%                 targets(row, :) = labels(d, :);
%             end
%         end
%         allOutputs = reshape(outputs', 1, numel(outputs));
%         targets = reshape(targets', 1, numel(allOutputs));
%         [tpr, fpr, thresholds] = roc(targets-1, allOutputs);
%         tpr = [tpr 1];
%         fpr = [fpr 1];   
%         plot(fpr, tpr, 'linewidth', 2, 'clipping', 'off', 'Color', colourSet(c,:));        
%         hold on
%         
%         for d=1:expSettings.nDatasets
%             datasetOutputs = combinedPostKnowns{d, i};
%             for r=1:nRepeats
%                 row = r+(d-1)*nRepeats;
%                 if sum(abs(outputs(row, :) - datasetOutputs(1, :, wmIdx))) > 0.5*expSettings.nSamples
%                     outputs(row, :) = 1-outputs(row, :);
%                 end
%             end
%         end
%         
%         allOutputs = reshape(outputs', 1, numel(outputs));
%         [tpr, fpr, thresholds] = roc(targets-1, allOutputs);
%         tpr = [tpr 1];
%         fpr = [fpr 1];    
%         plot(fpr, tpr, 'linewidth', 2, 'clipping', 'off', 'Color', colourSet(d, :));
%         hold all;
% 
%         axis([0 1 0 1]);
%         if d==1
%             legend(combMethods);
%         end
%         xlabel('False Positive Rate');
%         ylabel('True Positive Rate');
%         nKnown = expSettings.propKnown(selectedI(i)) * expSettings.nSamples;
%         title(sprintf('%s with %d known targets', combMethods{c}, nKnown));        
%     end
    set(h, 'PaperPositionMode', 'Auto');
    saveas(h, sprintf('%s/separate_%d', resultsDir, i), 'png');
    saveas(h, sprintf('%s/separate_%d', resultsDir, i), 'fig');
end

% ****

% 
% h = figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)*0.9]);
% ibccIdx = [1 2 3 10];
% %Just the IBCC variants
% for i=1:length(selectedI)
%     
%     p = selectedI(i);
%         for c=1:length(ibccIdx)
%             
%             if ibccIdx(c) > size(combinedPostKnowns{d,selectedI(i)}, 3)
%                 continue;
%             end
%             
%             subplot(length(selectedI), length(ibccIdx), c+(i-1)*length(ibccIdx));
%             for d=1:expSettings.nDatasets
%                 
%                 datasetOutputs = combinedPostKnowns{d, selectedI(i)};
%                 
%                 targets = zeros(1, size(datasetOutputs, 2));
%                 targets(1, :) = labels(d, (1:size(datasetOutputs, 2)));
%                 
%                 knownLabels = ExpRunner.getKnownTargets(...
%                     targets, ...
%                     expSettings.propKnown(p) * length(targets),...
%                     false);
%                 targets = targets(find(knownLabels==noScore));
%                 datasetOutputs = datasetOutputs(:, find(knownLabels==noScore), :);
%                 outputs = zeros(expSettings.nRepeats, size(datasetOutputs, 2));
%                 
%                 nRepeats = expSettings.nRepeats;
%                 if ibccIdx(c)==10
%                     for r=1:nRepeats
%                         post = datasetOutputs(r, :, ibccIdx(c));
%                         row = r;
%                                                 
%                         outputs(row, :) = post;
%                         [tpr, fpr, thresholds] = roc(targets-1, outputs(row,:));
%                         tpr = [tpr 1];
%                         fpr = [fpr 1];    
%                         plot(fpr, tpr, 'linewidth', 2, 'clipping', 'off', 'Color', colourSet(mod(d,size(colourSet,1))+1, :));                        
%                     end
%                     outputs = reshape(outputs', 1, numel(outputs));
%                     targets = targets' * ones(1, nRepeats);
%                     targets = reshape(targets, 1, numel(outputs));
%                     
%                 else
%                     post = datasetOutputs(1, :, ibccIdx(c));
%                     row = (d-1)*expSettings.nDatasets+1;
%                     outputs = post;
%                     [tpr, fpr, thresholds] = roc(targets-1, outputs);
%                     tpr = [tpr 1];
%                     fpr = [fpr 1];                    
%                     plot(fpr, tpr, 'linewidth', 2, 'clipping', 'off', 'Color', ...
%                         colourSet(mod(d, nColours)+1, :) );
%                 end
% 
%                 hold all;
%             end
%             axis([0 1 0 1]);
%             if d==1
%                 legend(combMethods);
%             end
%             xlabel('False Positive Rate');
%             ylabel('True Positive Rate');
%             nKnown = expSettings.propKnown(p) * expSettings.nSamples;
%             title(sprintf('%s with %d known targets', combMethods{ibccIdx(c)}, nKnown));        
%         end
% end
% 
% set(h, 'PaperPositionMode', 'Auto');
% saveas(h, sprintf('%s/ibccvariants_%d', resultsDir, i), 'png');
% saveas(h, sprintf('%s/ibccvariants_%d', resultsDir, i), 'fig');

% 
% %figure;
% % for i=1:length(selectedI)
% %     subplot(2, length(selectedI), i);
% %     for c=1:length(combMethods)       
% %        outputs = zeros(expSettings.nRepeats*expSettings.nDatasets, expSettings.nSamples);
% %        targets = zeros(expSettings.nRepeats*expSettings.nDatasets, expSettings.nSamples);
% %        
% %         for d=1:expSettings.nDatasets
% %             
% %             datasetOutputs = combinedPostKnowns{d, i};
% % 
% %             nRepeats = expSettings.nRepeats;
% %             for r=1:nRepeats
% %                 post = datasetOutputs(r, :, c);
% %                 row = r + (d-1)*nRepeats;
% %                 outputs(row, :) = post;   
% %                 targets(row, :) = labels(d, :);
% %             end
% %         end
% %         outputs = reshape(outputs', 1, numel(outputs));
% %         targets = reshape(targets', 1, numel(outputs));
% % 
% %         [tpr, fpr, thresholds] = roc(targets-1, outputs);
% %         tpr = [tpr 1];
% %         fpr = [fpr 1];    
% %         plot(fpr, tpr, 'linewidth', 2, 'clipping', 'off', 'Color', colourSet(c,:));        
% %         
% %         hold all;
% %         axis([0 1 0 1]);
% % 
% %         xlabel('False Positive Rate');
% %         ylabel('True Positive Rate');
% %         nKnown = expSettings.propKnown(selectedI(i)) * expSettings.nSamples;
% %         title(sprintf('ROC with %d known targets', nKnown));        
% %     end    
% % end
h = figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)/2.4]);
set(gca, 'FontSize', 16);
set(gca, 'FontName', 'Times New Roman');

auc = zeros(length(selectedI), length(combMethods));

for i=1:length(selectedI)
    subplot(1, length(selectedI), i);
        
    for c=1:length(combMethods)  
        
        p = selectedI(i);
        if exist('dataLengths', 'var')
            iPropKnown = expSettings.iPropKnown;
        else
            iPropKnown = selectedI(i);
        end
        nKnown = expSettings.propKnown(iPropKnown) * size(combinedPostKnowns{1, p}, 2);      
        nUnknown = (1-expSettings.propKnown(iPropKnown)) * size(combinedPostKnowns{1, p}, 2);
        
        outputs = zeros(expSettings.nRepeats*expSettings.nDatasets, nUnknown);
        targets = zeros(expSettings.nRepeats*expSettings.nDatasets, nUnknown);
       
        for d=1:expSettings.nDatasets
            
            datasetOutputs = combinedPostKnowns{d, p};

            nRepeats = expSettings.nRepeats;
            for r=1:nRepeats
                post = datasetOutputs(r, :, c);
                row = r + (d-1)*nRepeats;
                dLabels = labels(d, (1:size(post, 2)));
                
                knownLabels = ExpRunner.getKnownTargets(...
                    dLabels, ...
                    nKnown,...
                    false);
                targets(row,:) = dLabels(find(knownLabels==noScore));
                outputs(row,:) = post( find(knownLabels==noScore));
            end
        end
        
        for d=1:expSettings.nDatasets
            datasetOutputs = combinedPostKnowns{d, p};
            
            dLabels = labels(d, (1:size(datasetOutputs, 2)));
            knownLabels = ExpRunner.getKnownTargets(...
                    dLabels, ...
                    nKnown,...
                    false);
            
            for r=1:nRepeats
                row = r+(d-1)*nRepeats;
                %if sum(abs(outputs(row, :) - datasetOutputs(1,
                %find(knownLabels==noScore), wmIdx))) >
                %0.5*(1-expSettings.propKnown(p))*expSettings.nSamples
                if sum(abs(outputs(row, :) - targets(row,:) + 1)) > 0.5*nUnknown
                    outputs(row, :) = 1-outputs(row, :); % flip 
                    display(sprintf('c:%s, r:%d, d:%d, i:%d', combMethods{c}, r, d, i));
                end
            end
        end  
        
        outputs = reshape(outputs', 1, numel(outputs));
        targets = reshape(targets', 1, numel(outputs));

        [tpr, fpr, thresholds] = roc(targets-1, outputs);
        tpr = [tpr 1];
        fpr = [fpr 1];  
        
        x = targets.*outputs + (1-targets).*(1-outputs);
        y = (1-targets).*outputs + targets .* (1-outputs);
        %auc(i, c) = sjrob.auc(x,y);
        %auc(i, c) = trapz(fpr, tpr);
        
        %coloured lines instead of different shapes
        %plot(fpr, tpr, 'linewidth', 2, 'clipping', 'off', 'Color',
        %colourSet(c,:));    
        
        ms = 6;
        mc = 'k';
        col = [];
        switch c
            case 1
                lineStyle = '--g';
            case 2
                lineStyle = '-+';
                col = [0 0.6 1];
                mc = col;
                ms = 9;
            case 3 
                lineStyle = '--ko';
            case 4
                lineStyle = '-.';
                col = [0.5 0.2 0.2];
            case 5
                lineStyle = '-.k';
            case 6
                lineStyle = ':b';
            case 7
                lineStyle = '-r';
        end
        n=20;
        
        s = linspace(1,n,length(fpr));        % sampling independent variable
        sn = [1 (2:n-1)+randn(1,n-2)/n n];  % add some noise to avoid overlap
        xrs = interp1(s,fpr,sn,'nearest');    % downsample to n datapoints
        yrs = interp1(s,tpr,sn,'nearest');        
        
        if ~isempty(col )
            plot(xrs, yrs, lineStyle, 'LineWidth', 2, 'MarkerFaceColor', mc, 'MarkerSize', ms, 'Color', col);        
        else
            plot(xrs, yrs, lineStyle, 'LineWidth', 2, 'MarkerFaceColor', mc, 'MarkerSize', ms);        
        end
        
        set(gca, 'FontSize', 16);
        %set(gca, 'FontName', 'Times New Roman');
        hold all;
        axis([-0.005 1 0 1.005]);
        set(gca, 'FontSize', 16);
        %set(gca, 'FontName', 'Times New Roman');            
        if i==length(selectedI) && c==length(combMethods) 
%             legend({'Best Classifier', 'VB-IBCC-Seq', 'VB-IBCC', 'EM-ICC', 'Gibbs-IBCC', 'Simple Maj.', 'Weighted Maj.', 'Weighted Sum'}, 'Location', 'SouthEast');
%             legend({'VB-IBCC', 'EM-ICC', 'Gibbs-IBCC', 'Simple Maj.', 'Weighted Maj.', 'Weighted Sum'}, 'Location', 'SouthEast');
            legend({'Mean', 'VB-IBCC', 'VB-DynIBCC', 'Weighted Sum', 'DLR'}, 'Location', 'SouthEast');
        end
        xlabel('False Positive Rate');
        
        if i==1
            ylabel('True Positive Rate');
        end
        title(sprintf('%d Known Targets', nKnown));        
    end    
end

set(h, 'PaperPositionMode', 'Auto');
saveas(h, sprintf('%s/allflipped_%d', resultsDir, i), 'png');
saveas(h, sprintf('%s/allflipped_%d', resultsDir, i), 'fig');

% h = figure('Position',[1 scrsz(4)/2 scrsz(3)/3 scrsz(4)/2.1]);
% for a=1:expSettings.nAgents
%     for d=1:expSettings.nDatasets    
%         outputs = basePost(a, :, d);
%         targets = labels(d,:);
% 
%         [tpr, fpr, thresholds] = roc(targets-1, outputs);
%         tpr = [tpr 1];
%         fpr = [fpr 1];    
%         plot(fpr, tpr, 'linewidth', 2, 'clipping', 'off');        
% 
%         hold all;
%     end
%     axis([0 1 0 1]);
% 
%     xlabel('False Positive Rate');
%     ylabel('True Positive Rate');
%     nKnown = expSettings.propKnown(selectedI(i)) * expSettings.nSamples;
%     title(sprintf('ROC for Base Classifiers'));        
% end
% 
% 
% set(h, 'PaperPositionMode', 'Auto');
% saveas(h, sprintf('%s/base%d', resultsDir, i), 'png');
% saveas(h, sprintf('%s/base%d', resultsDir, i), 'fig');

% h = figure('Position',[1 scrsz(4)/2 scrsz(3)*1.2 scrsz(4)/2]);
% 
% selectedMethods = [1 2 3 8 9 10];
% 
% for i=1:length(selectedI)
%     subplot(1, length(selectedI), i);
%     for c=selectedMethods  
%         
%        outputs = zeros(expSettings.nRepeats*expSettings.nDatasets, size(combinedPostKnowns{1, i}, 2));
%        targets = zeros(expSettings.nRepeats*expSettings.nDatasets, size(combinedPostKnowns{1, i}, 2));
%        
%         for d=1:expSettings.nDatasets
%             
%             datasetOutputs = combinedPostKnowns{d, selectedI(i)};
% 
%             nRepeats = expSettings.nRepeats;
%             for r=1:nRepeats
%                 post = datasetOutputs(r, :, c);
%                 row = r + (d-1)*nRepeats;
%                 outputs(row, :) = post;   
%                 targets(row, :) = labels(d, (1:size(post, 2)));
%             end
%         end
%         
%         for d=1:expSettings.nDatasets
%             datasetOutputs = combinedPostKnowns{d, i};
%             for r=1:nRepeats
%                 row = r+(d-1)*nRepeats;
%                 if sum(abs(outputs(row, :) - datasetOutputs(1, :, wmIdx))) > 0.5*expSettings.nSamples
%                     outputs(row, :) = 1-outputs(row, :);
%                 end
%             end
%         end  
%         
%         outputs = reshape(outputs', 1, numel(outputs));
%         targets = reshape(targets', 1, numel(outputs));
% 
%         [tpr, fpr, thresholds] = roc(targets-1, outputs);
%         tpr = [tpr 1];
%         fpr = [fpr 1];    
%         plot(fpr, tpr, 'linewidth', 2, 'clipping', 'off', 'Color', colourSet(c,:));        
%         
%         hold all;
%         axis([0 1 0 1]);
%         if i==1
%             legendStrings = [combMethods(1:3) combMethods(8:10)];
%             legend(legendStrings, 'Location', 'SouthEast');
%         end
%         xlabel('False Positive Rate');
%         ylabel('True Positive Rate');
%         nKnown = expSettings.propKnown(selectedI(i)) * expSettings.nSamples;
%         title(sprintf('ROC with %d known targets, inverted results where mean error > 0.5', nKnown));        
%     end    
% end
% 
% 
% set(h, 'PaperPositionMode', 'Auto');
% saveas(h, sprintf('%s/mainflipped_%d', resultsDir, i), 'png');
% saveas(h, sprintf('%s/mainflipped_%d', resultsDir, i), 'fig');