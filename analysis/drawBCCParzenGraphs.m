%clear

if ~exist('threeDim', 'var')
    threeDim = true;
end

import settings.*;

vbIbccPaper.exp2seq;
%bccSettingsSimple1;
%bccLambdaSettings;
%bccSettings;
%bccSettingsHard;
%bccSettingsExp3;
%bccSettingsInversions;
%bccSettingsDynamicsExp3;

%uncomment this to draw the graphs for the lambda variation tests
%expSettings.multCombTestFile = sprintf('lambda%s', expSettings.multCombTestFile);

testDataFile = sprintf('%s%s_test_data.mat', expSettings.getDataDir(), ...
    expSettings.expLabel);
labelsFile = sprintf('%s%s_labels.mat', expSettings.getDataDir(), ...
    expSettings.expLabel);

if exist(testDataFile, 'file')
    labelledTestData = dlmread(testDataFile);
    labels = labelledTestData(:, expSettings.nSensors()+1);
elseif exist(labelsFile, 'file')
    labels = dlmread(labelsFile);
else
    display('cannot find labels');
end    

labels = reshape(labels, expSettings.nSamples, expSettings.nDatasets)';

combMethodsFile = sprintf('%scombMethods-%s', expSettings.getCombinerDir(), expSettings.multCombTestFile);
if exist(combMethodsFile, 'file')
    load(combMethodsFile);
else
    combMethods = {
        combiners.bcc.Ibcc.shortLabel ...
        }; 
end

%REMOVE THIS FOR NON_UPDATED RESULTS
%expSettings.multCombTestFile = sprintf('updatedNB_%s', expSettings.multCombTestFile);
load(sprintf('%s%s', expSettings.getCombinerDir(), expSettings.multCombTestFile));

if exist('combinedPostKnowns', 'var');
    varSettings = expSettings.propKnown .* expSettings.nSamples;

    %3d graphs
    if threeDim
        for c=1:length(combMethods)
            combLabel = combMethods{c};
                
            graphMaker = graphs.VariableSettingGraph(...
                sprintf('Errors of method "%s" with Varying no. Known Targets ', combLabel), ...
                expSettings.nSamples);    
            graphMaker.xLimit = [0 0.45*expSettings.nSamples];
         
            
            graphMaker.drawParzenGraph(true, combinedPostKnowns, labels, ...
                expSettings.nRepeats, ...
                varSettings, 0.01, 'No. Known Targets', c);
        end
    end
    
    %2d graph with all combination methods
    %probably want to collect data with low number of different settings
    nToPlot = 3; %length(expSettings.propKnown)
    plotIdx = 1;
    S = [1 2 6];
    for s=S
        graphMaker = graphs.VariableSettingGraph(...
            sprintf('Errors with %2.0f%% Targets Known', expSettings.propKnown(s)*100), ...
            expSettings.nSamples, length(combMethods));  
        graphMaker.fontsize = 24;
        for c=1:length(combMethods)
            
            if c==1 && s==1
                reuseGraph = false;
            else
                reuseGraph = true;
            end
            
            %ignore some unwanted methods - hack for report graphs only
%             if c == 9 || c == 8
%                 continue
%             end
%             
%             if s==1 &&(c==5 || c ==6)
%                 continue;
%             end
            
%             if c > 1 && length(combMethods) > 7
%                 names = [combMethods(1:3) combMethods(5:6) combMethods(8)];
%             elseif c > 1
%                 names = combMethods;
%             else
%                 names = [];
%             end
            
            if s<nToPlot
                graphMaker.labelXAxis = false;
            else
                graphMaker.labelXAxis = true;
            end
            
            if s==S(1)
                names = combMethods;%{'mean' 'maj. voting' 'DLR' 'N.B. LogOp'  'weighted sum' 'rounded ws' 'weighted maj.' 'unrounded wm' 'Ibcc'};
                graphMaker.dashed = false;
                if c==1
                    figure;
                end
            else
                %names = {'mean/weighed sum' 'majority voting/weighted majority' 'DLR' 'Ibcc'};
                names = [];
                graphMaker.dashed = true;
                %graphMaker.yLimit = [0 30];

            end
            
            subplot(nToPlot, 1, plotIdx);
            
            graphMaker.drawParzenGraph(true, combinedPostKnowns(:, s), labels, ...
                expSettings.nRepeats, varSettings(s), 0.01, ...
                sprintf('Errors with %2.2f%% Targets Known', expSettings.propKnown(s)*100), ...
                c, true,...
                names);
        end
        plotIdx = plotIdx + 1;
    end
    
elseif exist('combinedPostLambda', 'var')
    %use 1 - lambdasym to get the bias towards correctness rather than
    %incorrectness of classifiers
    varSettings = 1 - expSettings.lambdaSym;

    graphMaker = graphs.VariableSettingGraph('Symmetry of lambda', ...
        expSettings.nSamples);    
    
    combinedPostSyms = graphMaker.results3dTo2d(combinedPostLambda, 3, ...
        expSettings.nRepeats, expSettings.nSamples);
    graphMaker.drawParzenGraph(false, combinedPostSyms, labels, ...
        expSettings.nRepeats*size(combinedPostLambda, 3), ...
        varSettings, 0.01, 'bias of lambda toward correct outputs from classifiers');
    
    varSettings = expSettings.lambdaMag;

    graphMaker = graphs.VariableSettingGraph('Magnitude of lambda', ...
        expSettings.nSamples);  
    
    combinedPostMags = graphMaker.results3dTo2d(combinedPostLambda, 2, ...
        expSettings.nRepeats, expSettings.nSamples);
        
    graphMaker.drawParzenGraph(false, combinedPostMags, labels, ...
        expSettings.nRepeats*size(combinedPostLambda, 2), ...
        varSettings, 0.01, 'magnitude of lambda');
else
    display('no data to process');
end