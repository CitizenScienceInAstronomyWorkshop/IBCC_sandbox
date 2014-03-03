import settings.*

%bccSettingsInversions;
vbIbccPaper.exp2;

testDataFile = sprintf('%s%s_test_data.mat', ...
        expSettings.getDataDir(), expSettings.expLabel);

%store the outputs of the base classifiers. Managing this at this level
%because experiment runner doesn't handle multiple datasets.
basePostFile = sprintf('%s%s_base_post.mat', ...
        expSettings.getDataDir(), expSettings.expLabel);

if ~exist(expSettings.getDataDir(), 'file') || ~exist(basePostFile, 'file')
    display(['could not find the data file' basePostFile]);
end

labelledTestData = dlmread(sprintf('%s%s_test_data.mat', ...
    expSettings.getDataDir(), expSettings.expLabel));

labels = labelledTestData( ...
    :, ...
    expSettings.nSensors()+1);

if exist(basePostFile, 'file')
    basePost = dlmread(basePostFile);
    basePost = reshape(basePost, size(basePost, 1), ...
        expSettings.nSamples, expSettings.nDatasets);        
else
    basePost = zeros(expSettings.nAgents, expSettings.nSamples, ...
    expSettings.nDatasets);
end

nAgents = size(basePost, 1);
nAgentInsts = nAgents * expSettings.nDatasets;
nSamples = size(basePost, 2);
nTotalSamples = nSamples * expSettings.nDatasets;

%basePost = permute(basePost, [1 3 2]);
basePost = reshape(basePost, nAgents, nTotalSamples);


allLabels = reshape(labels', expSettings.nDatasets, nSamples);
labels = (labels * ones(1, nAgents))';
%labels = reshape(labels', [nAgents nSamples expSettings.nDatasets]);
%labels = permute(labels, [1 3 2]);
%labels = reshape(labels, [nAgentInsts nSamples]);

errors = abs(basePost - labels);
baseErrors = reshape(errors, nAgents*nSamples, expSettings.nDatasets);
meanErrors = zeros(nAgents, expSettings.nDatasets);
errorVars = zeros(nAgents, expSettings.nDatasets);

for d=1:expSettings.nDatasets
    meanErrors(:, d) = sum(errors(:, (d-1)*nSamples+1:d*nSamples), 2) ./ nSamples;
    for a=1:nAgents
        errorVars(a, d) = sum((errors(a, (d-1)*nSamples+1:d*nSamples) - meanErrors(a,d)).^2, 2) ./ nSamples;
    end
    for a=1:nAgents
        baseErrors(nSamples*(a-1)+1:nSamples*a, 1) = errors(a, (d-1)*nSamples+1:d*nSamples)';
    end
end

colourSet = graphs.createColourSet(nAgents);

% figure;
% for a=1:nAgents
%     plot(1:expSettings.nDatasets, meanErrors(a, :), 'Color', colourSet(a, :));
%     hold all;
% end
% hold off;

graphMaker = graphs.VariableSettingGraph('', nSamples*nAgents, expSettings.nDatasets);
graphMaker.fontsize = 24;
graphMaker.errorWindowSize = 1;
graphMaker.labelYAxis = false;
%graphMaker.yLimit = [0 2000];
graphMaker.xLimit = [0 1];
graphMaker.shade = true;
%The following produces a two-D graph (hard to read with 10 datasets)
figure;
legendStrings = [];

columns = 2;

for d=1:expSettings.nDatasets
    subplot(ceil(expSettings.nDatasets/columns), columns, d + mod(columns-expSettings.nDatasets, columns));
    
    if d > expSettings.nDatasets - columns
        graphMaker.labelXAxis = true;
    else
        graphMaker.labelXAxis = false;
    end
    
    legendStrings = [legendStrings num2str(d)];
    graphMaker.label = sprintf('Dataset %i', d);
    graphMaker.drawParzenGivenError(meanErrors(:, d)', errorVars(:, d)', 0.2, ...
        true, d, num2str(d), d);
    set(gca, 'XTick', [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
end

% graphMaker = graphs.VariableSettingGraph('Mean Errors of Base Classifiers', nSamples*nAgents, expSettings.nDatasets);
% 
% %The following produces a 3-D graph
% graphMaker.drawParzenGivenError(meanErrors', errorVars', 0.05, ...
%     false, 1:expSettings.nDatasets, 'Dataset');

% %The following produces a 3-D graph with perpendicular edges between
% %datasets
% varSettings = zeros(1, expSettings.nDatasets * 2);
% meanErrorsExp = [];
% errorVarsExp = [];
% for d=1:expSettings.nDatasets
%     meanErrorsExp = [meanErrorsExp meanErrors(:, d) meanErrors(:, d)];
%     errorVarsExp = [errorVarsExp errorVars(:, d) errorVars(:, d)];
%     varSettings(d*2 - 1) = d;
%     varSettings(d*2) = d + 1;
% end
% graphMaker.drawParzenGivenError(meanErrorsExp', errorVarsExp', 0.1, ...
%     false, varSettings, 'Dataset');
% 
% %draws the errors for all
% meanErrorsAll = reshape(meanErrors, nAgents*expSettings.nDatasets, 1);
% errorVarsAll = reshape(errorVars, nAgents*expSettings.nDatasets, 1);
% varSettings = 1;
% graphMaker.drawParzenGivenError(meanErrorsAll', errorVarsAll', 0.1, ...
%     false, varSettings, 'Dataset', 1);
% 
% %figure;
% %plot(reshape(ones(nAgents,1) * (1:expSettings.nDatasets), [1, expSettings.nDatasets*nAgents]),...
% %    reshape(meanErrors, [1, expSettings.nDatasets*nAgents]), '.');
% 
% agentsFile = sprintf('%s%s_multTestAgents.mat', ...
%         expSettings.getDataDir(), expSettings.expLabel);
% if exist(agentsFile, 'file')
%     load(agentsFile);
% 
%     informed = zeros(1, expSettings.nDatasets);
%     uninformed = zeros(1, expSettings.nDatasets);
% 
%     for a=1:expSettings.nAgents
%         for d=1:expSettings.nDatasets
%             if sum(agents{d}{a}.s(1:expSettings.nInfSensors)) > 0
%                 informed(d) = informed(d) + 1;
%             else
%                 uninformed(d) = uninformed(d) + 1;
%             end
%         end
%     end
% else
%     informed = sum(meanErrors<0.33, 1);
%     uninformed = sum(meanErrors >= 0.33, 1);
% end
% 
% figure;
% bar([informed' uninformed'], 'stacked');
% %plot(1:expSettings.nDatasets, informed);
% %hold all;
% %plot(1:expSettings.nDatasets, uninformed);
% legend('informed agents', 'uninformed agents');
% title('Number of Informed and Uninformed Agents in Each Dataset');
% hold off;
% 
% %graphMaker.drawParzenGivenError(baseErrors', 0.004, ...
% %    false, 1:expSettings.nDatasets, 'Dataset');