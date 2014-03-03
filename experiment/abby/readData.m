dataFile = '/homes/49/edwin/data/abby/extrafeats.dat';
% dataFile = '/Users/ablev/Documents/MATLAB/ibcc/code/experiment/solid/all/extrafeats.dat';

X = dlmread(dataFile, '\t', 1, 0); % order is: source_id  test_id prob(0) prob(1) gold_lbl

numSources = X(end,1); % aka nAgents
numTotPts = X(end,2);
numTstPts = 100;
numClasses = 2;

% get 1-indexed gold labels
goldlbls = X(1:numTotPts,5)';

[rows cols] = size(goldlbls);
trainlbls = goldlbls;
trainlbls(:, cols-numTstPts+1:cols) = 0;
cellxx ={X(:,1), X(:,2), round(X(:,4))}; % used in plain and mixed IBCC 
sp = sparse(X(:,1), X(:,2), X(:,4), numSources, numTotPts); % used in CIBCC

bccSettings = settings.BccSettings();
% singleAlpha = [1 1; 1 1]; 
singleAlpha = [20 4; 4 20]; 
bccSettings.nu = {[100 100]};
bccSettings.Alpha = repmat(singleAlpha, [1,1,numSources]);

%%%%%%% Plain IBCC %%%%%%
%disp('IBCC');
%combiner = combiners.bcc.IbccVb(bccSettings, numSources, 1, trainlbls, [], numClasses, numClasses);
%[combinedPost Alpha] = combiner.combineDecisions(cell);
%disp('P(1) of test points:');
%combinedPost(:, (end-numTstPts+1):end)
%goldlbls(:, (end-numTstPts+1):end)

%%%%%%% Conintuous IBCC %%%%%%
disp('CIBCC:');
combiner = combiners.bcc.CIbccVb(bccSettings, numSources, 1, trainlbls, [], numClasses, numClasses);
combiner.useLikelihood = false;
[combinedPost Alpha] = combiner.combineDecisions(sp);
disp('P(1) of test points:');
testResults = combinedPost(:, (end-numTstPts+1):end)
testlbls = goldlbls(:, (end-numTstPts+1):end)

[auc,fig,thresholds] = graphs.ClassifierPerformanceGraph.drawRoc(testResults, testlbls-1, ...
    {'CIBCC'}, false, false, true);
auc

[auc,fig,thresholds] = graphs.ClassifierPerformanceGraph.drawRoc(round(testResults), testlbls-1, ...
    {'CIBCC'}, false, false, true);
auc

% return; 

%%%%%%% Mixed IBCC (dynamic + continuous IBCC) %%%%%%
disp('MIBCC');
%MixedIbccVb(bccSettings, nFeat, featSettings, features, nAgents, K, targets, agents, nClasses, nScores)
combiner = combiners.bcc.MixedIbccVb(bccSettings, numSources, bccSettings, sp, numSources, 1, trainlbls, [], numClasses, numClasses);
combiner.initAll = true;
[combinedPost Alpha] = combiner.combineDecisions(cellxx);
disp('P(1) of test points:');
combinedPost(:, (end-numTstPts+1):end)
goldlbls(:, (end-numTstPts+1):end)

[auc,fig,thresholds] = graphs.ClassifierPerformanceGraph.drawRoc(testResults, testlbls-1, ...
    {'CIBCC'}, false, false, true);
auc

[auc,fig,thresholds] = graphs.ClassifierPerformanceGraph.drawRoc(round(testResults), testlbls-1, ...
    {'CIBCC'}, false, false, true);
auc

%{ 
xnew = [,1:numPredPts]';
for i = 1:numSources
    idx = i * numPredPts;
    tmp = X(X(:,1)==i,:);
    xnew = [xnew tmp(:,4)];
end
%}
