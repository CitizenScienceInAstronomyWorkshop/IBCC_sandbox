%  clear
%use different settings files for each experiment.
% rmdir('/homes/49/edwin/data/thesis_synth/exp1');
settingsFuncs = { ...
...%     @settings.thesis.exp1, ...
    ...%@settings.thesis.exp2, ...
    ...%@settings.thesis.exp3, ...
    ...%@settings.thesis.exp7, ...
    @settings.thesis.exp4, ...
    ...%@settings.thesis.exp6, ...
    };

display([ 'Number of experiments: ' num2str(numel(settingsFuncs)) ]);

for expNo = 1:numel(settingsFuncs)
    
    varySettingFunc = settingsFuncs{expNo};

    [expSet, synthSet, bccSet, ~, range, varName] = varySettingFunc(1);
    outputDir = expSet.outputDir;
    
    nSettings = length(range);

    if ~exist(expSet.getCombinerDir(), 'file')
        mkdir(expSet.getCombinerDir());
    end

    %pick the combination methods to use.
    combMethods = {...
...%         combiners.MeanDecision.shortLabel,... %1 - blue
...%           combiners.bcc.IbccVb.shortLabel,... %2 - blue
...%            combiners.bcc.IbccSampling.shortLabel, ... %2 - green
...%           combiners.DlrCombiner.shortLabel, ... % 1 - green
...%            combiners.weighted.WeightedSum.shortLabel,... % 1 - red
           combiners.SimpleMajorityVoting.shortLabel,... %2 - blue
...%            combiners.weighted.WeightedMajority.subShortLabel, ... %2 s
...%           combiners.weighted.NaiveBayesLogodds.shortLabel,... %1
          %combiners.weighted.Bma.subSubShortLabel,... %1
         %use these in comparison of IBCC inference methods in next chapter
    %     combiners.bcc.IbccVb.shortLabel,...
        };
    combLabels = {...
...%         'Mean',...
...%         'IBCC-VB', ...
         'IBCC-Gibbs',...
...%         'DLR',...
         'Weighted Sum',...
        'Maj. Voting',...
        'Weighted Maj.',...
...%         'NB-LogOP',...
        };
    
    graph1 = [1 2 3 4];%[3 4 5 8];
    graph2 = [];%[1 2 3 6 7];
    
    graphColors = {[0 0.8 0.8], [0 0 1], [0 0.6 0], [0.5 0.2 0.8], [1 0 0], [0.8 0.8 0], [1 0 1], [1 0.5 0], [0 0 0]};
    graphMarkers = {'-','-s','-*','-s','-.o','-x','->','-d','--+'};
    
    graph1Colors = graphColors([graph1 end]);
    graph1Markers = graphMarkers([graph1 end]);
    graph2Colors = graphColors([graph2 end]);
    graph2Markers = graphMarkers([graph2 end]);

    rerunAgentsEachIt = true;
    multCombTests;
    
    %Show the following in thesis:

    %1. Graphs for each method separately, showing their changes with the
    %independent variable.
    %2. Comparative graph of 4 important methods (weighted sum, IBCC, NB logop,
    %simple maj) for the middle variable setting (just an extra visual summary
    %of the previous graphs
    %3. Table of AUCs. 
    %4. Graph of AUCs against the independent variable. 
    % This totals 11 graphs and one table per experiment. Alternatively, show
    % comparative graphs for each variable setting instead of (1) and (2) and
    % get 5+1=6 graphs in total for each experiment.

    %Remember the IBCC VB lines don't belong here. 
    % There may be some experiments where some of the methods do not need to be
    % compared and graphs can be simplified.

    if strcmp(varName, 'Mean Error Rate for Class 1')
        range = errorSum;
    end
        
    
    %show Auc
    close all
    if exist('noGraphsPlease','var')
        clear progressMonitor
        clear combinedPostKnowns
        clear labelsKnowns
        continue
    end
    
    [aucs,meanAucs,devAucs,agentAucs] = drawThesisGraphs(synthSet.nDatasets, expSet, range, ...
            varName, combinedPostKnowns, agentPostKnowns, labelsKnowns, combLabels, nSettings, ...
            graphColors, graphMarkers, graph1, graph2, graph1Colors, graph1Markers, graph2Colors, graph2Markers);
    
    
%     if expNo==1
%         close all
%         %redraw graphs with agent error as the x axis
%         range = aucs(:,end);
%         aucs = drawThesisGraphs(synthSet.nDatasets, expSet, range, ...
%             varName, combinedPostKnowns, agentPostKnowns, labelsKnowns, combMethods, nSettings);        
%     end
%     
    %show Bayes error 
    [errors] = drawThesisError('Brier Score (Mean Squared Error)', expSet, ...
        range, varName, combinedError, bestAgentError, meanAgentError, combLabels, ...
        graph1, graph2, graph1Colors, graph1Markers, graph2Colors, graph2Markers);
%     fprErrors = drawThesisError('False Positive Rate', expSet, range, varName, fpr, bestAgentFpr, combLabels, graph1, graph2);
%     fnrErrors = drawThesisError('False Negative Rate', expSet, range, varName, fnr, bestAgentFnr, combLabels, graph1, graph2);
%     entropies = drawThesisError('Entropy Over Target Labels', expSet, range, varName, entropy, [], combLabels, graph1, graph2);
    
    %save data for tables   
    dlmwrite([outputDir '/experiment_' num2str(expNo) '_meanaucs.csv'], meanAucs);
    dlmwrite([outputDir '/experiment_' num2str(expNo) '_devaucs.csv'], devAucs);
    dlmwrite([outputDir '/experiment_' num2str(expNo) '_meanerrors.csv'], errors);
    
    devErrors = sum((repmat(combinedError, [1,1,synthSet.nDatasets])-combinedErrorAll).^2, 3) ./ synthSet.nDatasets;
    devErrors = devErrors .^ 0.5;
    
    dlmwrite([outputDir '/experiment_' num2str(expNo) '_deverrors.csv'], devErrors);

    save([outputDir '/experiment_' num2str(expNo) '_progressmon.mat'], 'progressMonitor');
    clear progressMonitor
    clear combinedPostKnowns
    clear labelsKnowns
end