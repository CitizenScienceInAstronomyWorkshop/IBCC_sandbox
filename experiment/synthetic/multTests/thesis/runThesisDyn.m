%  clear
%use different settings files for each experiment.
% rmdir('/homes/49/edwin/data/thesis_synth/exp1');
settingsFuncs = { ...
         @settings.thesis.dyn11,...
%         @settings.thesis.dyn5,...
%         @settings.thesis.dyn9,...
...%         @settings.thesis.dyn3,...
...%         @settings.thesis.dyn2,...
...%         @settings.thesis.dyn1,...
...%         @settings.thesis.dyn4,...
...%         @settings.thesis.dyn6,...
...%             @settings.thesis.dyn7,...
...%             @settings.thesis.dyn8,...
...%             @settings.thesis.dyn10,...
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
...        combiners.MeanDecision.shortLabel,... %1 - blue
        combiners.bcc.IbccVb.shortLabel,... %2 - blue
        combiners.bcc.DynIbccVb.shortLabelInitAll,...
        combiners.bcc.IbccSampling.shortLabel,...
...%        combiners.DlrCombiner.shortLabel, ... % 1 - green
...%        combiners.weighted.WeightedSum.shortLabel,... % 1 - red
...%        combiners.weighted.WeightedMajority.subShortLabel, ... %2 s
        };

    display('remember we can also change the initialisation of DynIBCC');

    combLabels = {...
        ...%'Mean',...
         'IBCC-VB', ...
         'DynIBCC-VB',...
         'IBCC-Gibbs',...
        ...%'DLR',...
        ...%'Weighted Sum',...
        ...%'Weighted Maj.',...
        };
    
    graph1 = [1 2 3];% 4 5 6];
    graph2 = [];
    
    graphColors = {[0 0.8 0.8], [0 0 1], [0 0.6 0], [0.5 0.2 0.8], [1 0 0], [0.8 0.8 0], [1 0 1], [1 0.5 0], [0 0 0]};
    graphMarkers = {'-<','-s','-*','-s','-.o','-x','->','-d','--+'};
    
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
    [aucs,meanAucs,devAucs,agentAucs] = drawThesisGraphs(synthSet.nDatasets, expSet, range, ...
        varName, combinedPostKnowns, agentPostKnowns, labelsKnowns, combLabels, nSettings, ...
        graphColors, graphMarkers, graph1, graph2, graph1Colors, graph1Markers, graph2Colors, graph2Markers);
    
    %show Brier score
    [errors] = drawThesisError('Brier Score (Mean Squared Error)', expSet, ...
        range, varName, combinedError, bestAgentError, meanAgentError, combLabels, ...
        graph1, graph2, graph1Colors, graph1Markers, graph2Colors, graph2Markers);
    
    %save data for tables   
    dlmwrite([outputDir '/experiment_' num2str(expNo) '_meanaucs.csv'], meanAucs);
    dlmwrite([outputDir '/experiment_' num2str(expNo) '_devaucs.csv'], devAucs);
    dlmwrite([outputDir '/experiment_' num2str(expNo) '_meanerrors.csv'], errors);
    
    devErrors = sum((repmat(combinedError, [1,1,synthSet.nDatasets])-combinedErrorAll).^2, 3) ./ synthSet.nDatasets;
    devErrors = devErrors .^ 0.5;
    
    dlmwrite([outputDir '/experiment_' num2str(expNo) '_deverrors.csv'], devErrors);

    save([outputDir '/experiment_' num2str(expNo) '_progressmon.mat'], 'progressMonitor');
    save([outputDir '/experiment_' num2str(expNo) '_synthSet.mat'], 'synthSet');
    clear progressMonitor
    clear combinedPostKnowns
    clear labelsKnowns
end