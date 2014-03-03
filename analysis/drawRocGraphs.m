%clear

if ~exist('threeDim', 'var')
    threeDim = true;
end

import settings.*;

bccSettings; expSettings.topSaveDir = '/homes/49/edwin/matlab/combination/data/vbPrsMerged/bcc3';
% bccSettingsHard; expSettings.topSaveDir = '/homes/49/edwin/matlab/combination/data/vbPrsMerged/bcc4Hard';
% bccSettingsExp3; expSettings.topSaveDir = '/homes/49/edwin/matlab/combination/data/vbPrsMerged/bccExp3';
% bccSettingsInversions; expSettings.topSaveDir = '/homes/49/edwin/matlab/combination/data/vbPrsMerged/bccInv';
% bccSettingsDynamicsExp3; expSettings.topSaveDir = '/homes/49/edwin/matlab/combination/data/vbPrsMerged/bccExp5';

[labels, combinedPostKnowns, combMethods] = loadMultTestResults(expSettings);

colourSet = graphs.createColourSet(length(combMethods));

colourSet(4, :) = [1 0.7 0];
%obj.colourSet(5, :) = [0.9 0.9 0];
colourSet(6, :) = [0.1 0.1 1];
colourSet(8, :) = [0 1 1];
colourSet(9, :) = [1 0.3 0.7];
colourSet(10, :) = [0.5 0 0.5];

%indexes of proportions of known targets that we want to use
selectedI = [1 2 6];
%draw a separate window for each number of known targets
for i=1:length(selectedI)
    figure;
        for d=1:expSettings.nDatasets
            
            subplot(2, ceil(expSettings.nDatasets/2), d);
            
            for c=1:length(combMethods)
                outputs = zeros(expSettings.nRepeats, expSettings.nSamples);
                targets = zeros(1, length(labels));
                targets(1, :) = labels(d, :);

                datasetOutputs = combinedPostKnowns{d, i};
            
                nRepeats = expSettings.nRepeats;
                if c==10
                    for r=1:nRepeats
                        post = datasetOutputs(r, :, c);
                        row = r;
                        outputs(row, :) = post;
                    end
                    outputs = reshape(outputs', 1, numel(outputs));
                    targets = targets' * ones(1, nRepeats);
                    targets = reshape(targets, 1, numel(outputs));
                    
                else
                    post = datasetOutputs(1, :, c);
                    row = (d-1)*expSettings.nDatasets+1;
                    outputs = post;
                end
                [tpr, fpr, thresholds] = roc(targets, outputs);
                tpr = [tpr 1];
                fpr = [fpr 1];

                plot(fpr, tpr, 'linewidth', 2, 'clipping', 'off', 'Color', colourSet(c, :));
                hold all;
            end
            axis([0 1 0 1]);
            if d==1
                legend(combMethods);
            end
            xlabel('False Positive Rate');
            ylabel('True Positive Rate');
            nKnown = expSettings.propKnown(selectedI(i)) * expSettings.nSamples;
            title(sprintf('%s with %d known targets', combMethods{c}, nKnown));        
        end

    %puts separate plots for each method
    %draw a separate subplot for each combination method
%     for c=1:length(combMethods)
%         subplot(2, ceil(length(combMethods)/2), c);
%         outputs = zeros(expSettings.nDatasets*expSettings.nRepeats, expSettings.nSamples);
%         for d=1:expSettings.nDatasets
%             targets = zeros(1, length(labels));
%             targets(1, :) = labels(d, :);
% 
%             datasetOutputs = combinedPostKnowns{d, i};
%             for r=1:expSettings.nRepeats
%                 post = datasetOutputs(r, :, c);
%                 row = (d-1)*expSettings.nDatasets+r;
%                 outputs(row, :) = post;
%                 
%                 [tpr, fpr, thresholds] = roc(targets, outputs(row, :));
%                 tpr = [tpr 1];
%                 fpr = [fpr 1];
%                 plot(fpr, tpr, 'linewidth', 2, 'clipping', 'off');
%                 hold all;
%             end
%         end
%         %axis([0 1 0 1]);
%         xlabel('False Positive Rate');
%         ylabel('True Positive Rate');
%         nKnown = expSettings.propKnown(selectedI(i)) * expSettings.nSamples;
%         title(sprintf('%s with %d known targets', combMethods{c}, nKnown));        
%     end
end