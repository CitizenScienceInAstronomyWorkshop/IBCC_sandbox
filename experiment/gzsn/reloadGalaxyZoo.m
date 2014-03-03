
% [snBaseOutputs, labels, typeLabels, classLabels, typeOnlyLabels, typeAssets, classAssets] = reloadGalaxyZooData(true, true);
% 
% import settings.*
% gz.gz_sn1;
% 
% %pick the combination methods to use.
% combMethods = {...
%     combiners.MeanDecision.shortLabel,...  
% ...%     combiners.SimpleMajorityVoting.shortLabel, ...
% ...%     combiners.WeightedMajority.subShortLabel, ...
% ...%     combiners.WeightedSum.shortLabel,...
%     combiners.IbccVb.shortLabel,...   
% ...%     combiners.IbccVb.shortLabelSeq, ...      
% %     combiners.Ibcc.shortLabel, ...
%     };
% 
% %put the labels from 7 and 8 into the main body of data

% ptfTypeAgent = max(snBaseOutputs{1}) + 1;
% ptfClassAgent = ptfTypeAgent + 1;
% nTypeKnown = length(typeLabels);
% nClassKnown = length(classLabels);
% 
% % Rebalance the class proportions 
% posAssets = find(labels==2);
% posAgents = sparse(max(snBaseOutputs{1}), 1);
% for p=posAssets'
%     agents = snBaseOutputs{1}(snBaseOutputs{2}==p);
%     posAgents(agents) = posAgents(agents) + 1;
% end
% 
% minAgentResp = 50;
% minFreqAgents = 10;
% 
% keyAgents = find(posAgents>=minAgentResp);
% 
% negAssets = find(labels==1);
% negAssetsToKeep = [];
% negAgents = sparse(max(snBaseOutputs{1}), 1);
% for n=negAssets'
%     display(n);
%     agents = unique(snBaseOutputs{1}(snBaseOutputs{2}==n));
%     agentCounts = posAgents(agents);
%     if sum(agentCounts>=minAgentResp) > minFreqAgents
%         negAssetsToKeep = [negAssetsToKeep; n];
%         negAgents(agents) = negAgents(agents) + 1;
%     end
% end
% 
% negPosImbalance = posAgents - negAgents;
% imbalanceIdxs = find(abs(negPosImbalance)>10);
% 
% assetsToKeep = [negAssetsToKeep; posAssets];
% 
% %subsample labels
% subLabels = sparse(length(labels), 1);
% subLabels(assetsToKeep) = labels(assetsToKeep);
% 
% %subsample decision data from base classifiers
% idxsToKeep = ismember(snBaseOutputs{2}, assetsToKeep);
% testData = {};
% testData{1} = snBaseOutputs{1}(idxsToKeep);
% testData{2} = snBaseOutputs{2}(idxsToKeep);
% testData{3} = snBaseOutputs{3}(idxsToKeep);
% %%%%%%%%%%%%%
% 
% nonZeroLabels = find(subLabels~=0);
% nLabels = numel(nonZeroLabels);
% 
% nTypeOnlyKnown = length(typeOnlyLabels);
% 
% typeAssets = double(typeAssets);
% classAssets = double(classAssets);
% 
% % % %run with k-fold cross validation
% % nFolds = 5; 
% % c_class = cvpartition(nLabels, 'kfold', nFolds);
% % c_type = cvpartition(nTypeKnown, 'kfold', nFolds);
% 
% nAgents = max(snBaseOutputs{1});
% nAssets = max(snBaseOutputs{2});

% for nRatio=[0]

    combinedPostAllFolds = [];
    labelsAllFolds = [];
    display(['nRatio: ' num2str(nRatio)]);

    %this is irrelevant now we're only looking at labelled data points
%     testData = selectPortionUnknown(nLabels, nRatio, snBaseOutputs, nonZeroLabels);  
    
    progressMonitor = scoring.CombinerProgressMonitor([], length(combMethods));
    progressMonitor.method = 'absolute error';
    
    totalTime = zeros(1, length(combMethods));    
    
    for i=1:nFolds
        
        display(['num folds: ' num2str(i)]);
    
        nClassTraining = sum(c_class.training(i));
%         nTypeOnlyTraining = sum(c_typeOnly.training(i));
    
        foldData = cell(length(testData), 1);
    
        testIdxs = nonZeroLabels(c_class.test(i));
        trainingIdxs = nonZeroLabels(c_class.training(i));      
        typeTrainIdxs = ~ismember(typeAssets, testIdxs);
        typeTrainAssets = typeAssets(typeTrainIdxs);        
        
        %data 1 is agents
        foldData{1} = [testData{1}; ... %data from human classifiers
%             ptfTypeAgent*ones(length(typeTrainIdxs), 1); ... %type data from ptf (assets with a ptf class)
            ];%ptfClassAgent*ones(nClassTraining, 1)]; %class data from ptf

        %data 2 is data point no.
        foldData{2} = [testData{2}; ...
%             typeTrainAssets; ...
            ];%classAssets(c_class.training(i))];

        %data 3 is scores given by agents
        foldData{3} = [testData{3}; ...
%             typeLabels(typeTrainIdxs); ...
            ];%classLabels(c_class.training(i))];

        display(['length of fold data=' num2str(length(foldData{1}))]);
                
        trainingLabels = subLabels;
        trainingLabels(testIdxs) = 0;
        testLabels = subLabels(testIdxs);
        
        expSettings.propKnown = 0;
        expSettings.iPropKnown = 1;

        %training labels are passed in to the runner - these are given to
        %combiners as known labels. The complete set of labels is then used
        %for evaluating results.
        runner = GzRunner(expSettings, foldData, trainingLabels, subLabels);
        progressMonitor.updateLabels(subLabels, testIdxs, length(trainingIdxs));
        runner.progressMonitor = progressMonitor;
       
        %new_data, runBaseAgents, drawGraphs, combMethods, sort results)
        combinedPost = runner.runCombineAll(false, 'dontTouch', false, combMethods, false); 
        
        totalTime = totalTime + runner.runTime;

        if exist('keepAllLabels', 'var') && keepAllLabels %...even the 
            testCombinedPost = combinedPost(:, trainingIdxs);
            testLabels = subLabels(trainingIdxs);
        else
            testCombinedPost = combinedPost(:, testIdxs);        
        end
        combinedPostAllFolds = [combinedPostAllFolds testCombinedPost];
        labelsAllFolds = [labelsAllFolds testLabels'];

        allFoldsDir = '/homes/49/edwin/matlab/combination/data/combination3_sandbox_';        
        
        if ~exist(sprintf('%sall%dFolds', allFoldsDir, nFolds), 'dir');
            mkdir(sprintf('%sall%dFolds', allFoldsDir, nFolds));
        end
                
        save(sprintf(...
            '%sall%dFolds/combin_minPosR%d_minFreqA%d_rUnk%1.1f.mat', ...
            allFoldsDir, nFolds, minAgentResp, minFreqAgents, nRatio), 'combinedPostAllFolds');    
        save(sprintf(...
            '%sall%dFolds/labels_minPosR%d_minFreqA%d_rUnk%1.1f.mat', ...
            allFoldsDir, nFolds, minAgentResp, minFreqAgents, nRatio), 'labelsAllFolds');    
        save(sprintf(...
            '%sall%dFolds/cvClass_minPosR%d_minFreqA%d_rUnk%1.1f.mat', ...
            allFoldsDir, nFolds, minAgentResp, minFreqAgents, nRatio), 'c_class');    
%         save(sprintf(...
%             '%sall%dFolds/cvType_minPosR%d_minFreqA%d_rUnk%1.1f.mat', ...
%             allFoldsDir, nFolds, minAgentResp, minFreqAgents, nRatio), 'c_type');    
        save(sprintf(...
            '%sall%dFolds/progMon_minPosR%d_minFreqA%d_rUnk%1.1f.mat', ...
            allFoldsDir, nFolds, minAgentResp, minFreqAgents, nRatio), 'progressMonitor');         
        
    end
    
    avgTime = totalTime ./ nFolds
    
%     
% end