
if ~exist('expSettings', 'var')
    import settings.*
    gz.gz_sn1;
end

%nRatio=[0 0.5 1 10]
nRatio = 0;

if ~exist('resultsAllFolds','var')
    load(sprintf('%sall%dFolds/combin_minPosR%d_minFreqA%d_rUnk%1.1f.mat', ...
            allFoldsDir, nFolds, minAgentResp, minFreqAgents, nRatio));    
end
if ~exist('labelsAllFolds', 'var')
    load(sprintf('%sall%dFolds/labels_minPosR%d_minFreqA%d_rUnk%1.1f.mat', ...
            allFoldsDir, nFolds, minAgentResp, minFreqAgents, nRatio));
end
if ~exist('nFolds', 'var')
    nFolds = 5;
end
if ~exist('c_class', 'var') && nFolds > 1
    load(sprintf('%sall%dFolds/cvClass_minPosR%d_minFreqA%d_rUnk%1.1f.mat', ...
            allFoldsDir, nFolds, minAgentResp, minFreqAgents, nRatio));
end

testOnlyResults = [];
testOnlyLabels = [];

%uncomment if the results contain all data points - this code selects only
%the test data, i.e. the data points we had no training label for but for
%which we do have a testing label
% for i=1:nFolds
%     
%     testIdxs = classAssets(c_class.test(i));
%     
%     testOnlyResults = [testOnlyResults resultsAllFolds(:, testIdxs)];
%     testOnlyLabels = [testOnlyLabels labelsAllFolds(:, testIdxs)];
% end
% 
% resultsAllFolds = testOnlyResults;
% labelsAllFolds = testOnlyLabels;

if ~exist('combMethods', 'var')
    %pick the combination methods to use.
    combMethods = {...
        combiners.MeanDecision.shortLabel,...
    ...%     combiners.bcc.IbccVb.shortLabelSeq, ...    
        combiners.SimpleMajorityVoting.shortLabel, ...
        combiners.weighted.WeightedMajority.subShortLabel, ...
        combiners.weighted.WeightedSum.shortLabel,...
        combiners.bcc.IbccVb.shortLabel,...   
        combiners.bcc.Ibcc.shortLabel, ...
        };
end

graphs.ClassifierPerformanceGraph.drawRoc(resultsAllFolds, labelsAllFolds-1, combMethods, false);

% legendStrings = cat(2, combMethods, 'real answers');
% graph = graphs.ClassifierPerformanceGraph(expSettings, {}, 1, legendStrings);
%                            
% combinedWithAnswers = [resultsAllFolds; labelsAllFolds-1];
%                 
% [sortedMeans, idxs] = sort(combinedWithAnswers(1,:));
% sortedLabels = labelsAllFolds(:, idxs);
%                 
% roundedResults = round(resultsAllFolds(2,:));
% errors = roundedResults - labelsAllFolds + 1;
% nFP = errors&roundedResults;
% nTP = (errors==0)&roundedResults;
% 
% % graph.drawResultHistogram(sortedLabels, 'Histogram of scores', {'not SN' 'SN'});
% graph.drawPostGraph(combinedWithAnswers(:, idxs), 'Output of Combination Methods');
% graph.drawErrorGraph(labelsAllFolds(:,idxs), resultsAllFolds(:,idxs), 'Absolute Error of Combination Methods');
%                 
% hfigs = get(0, 'children');  %Get list of figures
% for g=1:length(hfigs)
%     filename = sprintf('%s/graphs/%s_graph%d', rootDir, datestr(now, 'yy_mm_dd__HH_MM_SS'), g);
%     saveas(hfigs(g), [filename '.fig']) %Matlab .FIG file
%     saveas(hfigs(g), [filename '.png'])
% end