%use different settings files for each experiment.
% rmdir('/homes/49/edwin/data/thesis_synth/exp1');

clear

settingsFuncs = {...
      @settings.thesis.exp1, ...
      @settings.thesis.exp2, ...
      @settings.thesis.exp3, ...
      @settings.thesis.exp4, ...
      @settings.thesis.exp6, ...
      @settings.thesis.exp7
    };

display([ 'Number of experiments: ' num2str(numel(settingsFuncs)) ]);

for expNo = 1:numel(settingsFuncs)
    
    varySettingFunc = settingsFuncs{expNo};

    [expSet, synthSet, bccSet, ~, range, varName] = varySettingFunc(1);
    outputDir = expSet.outputDir;
    rootSuffix = 'VB/';
    replacementDir = [expSet.rootDir rootSuffix];
    expSet.setDirNames(replacementDir);
    
    nSettings = length(range);

    if ~exist(expSet.getCombinerDir(), 'file')
        mkdir(expSet.getCombinerDir());
    end

    %pick the combination methods to use.
    combMethods = {...
...%         combinerrs.bcc.IbccSampling.shortLabel, ... %2
        combiners.bcc.IbccEm.shortLabel,... %2        
        combiners.bcc.IbccVb.shortLabel,...
        };
    
    combLabels = {...
        'IBCC-Gibbs', ... %2
        'IBCC-VB',...
        };
    
    specialExp = 3;
    graph1 = [2 1];
    graph2 = [];

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
    aucs = drawThesisGraphs(synthSet.nDatasets, expSet, range, ...
        varName, combinedPostKnowns, agentPostKnowns, labelsKnowns, combLabels, nSettings, graph1, graph2);
    
%     if expNo==1
%         close all
%         %redraw graphs with agent error as the x axis
%         range = aucs(:,end);
%         aucs = drawThesisGraphs(synthSet.nDatasets, expSet, range, ...
%             varName, combinedPostKnowns, agentPostKnowns, labelsKnowns, combMethods, nSettings);        
%     end
%     
    %show Bayes error 
    errors = drawThesisError('Brier Score (Mean Squared Error)', expSet, range, varName, combinedError, errorSum, combLabels, graph1, graph2);
    
    %save data for tables   
    dlmwrite([outputDir '/experiment_' num2str(expNo) '_aucs.csv'], aucs);
    dlmwrite([outputDir '/experiment_' num2str(expNo) '_errors.csv'], errors);
end