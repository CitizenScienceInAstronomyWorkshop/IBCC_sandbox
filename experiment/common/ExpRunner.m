classdef ExpRunner < handle
    
    properties        
        labelledTestData %for the base classifiers, not for the combiners
        labelledTrainingData % for the base classifiers, not the combiners
        labels %complete set of labels
        trainingLabels = [];%for the combiners
        
        agents
        basePost = [];
        baseLogodds = [];
        
        combinedPost = [];
        %as an alternative to the basePost, some combiners may take an Nx3
        %matrix with column 1 containing agent IDs, column 2 containing
        %data-point IDs, and column 3 containing the score/posterior
        %produced by a base agent. This allows the same agents to provide 
        %scores to the same object more than once. 
        % Exp should convert between the basePost and baseScoreSet
        % forms for each classifier.
        baseScoreSet = [];
        
        agentConfidence
        
        changes
        
        expSet
        bccSet
        synthSet
        
        clusterer 
        
        combiner        
        
        loadData = true;
        saveData = false;
        
        %runBaseClassifiers = true;

        convertTargetsToAgents = false;
        
        progressMonitor = [];
        
        runTime = [];
    end    
    
    methods (Static)
        function targets = getKnownTargets(labels, nKnown, spreadTargets)
            targets = sparse(ones(1,length(labels)), 1:length(labels), labels);

            if spreadTargets 
                
                %spread randomly
                display('spreading training labels randomly');
                idxs = randsample(length(labels), length(labels)-nKnown); %idxs for the UNKNOWN 
                targets(idxs) = 0;
                
                % spread the known labels into 5 lumps
%                 nBurst = round(nKnown / 5);
%                 nPeriod = round(length(labels) / 5);
% 
%                 for n=1:length(labels)
%                     if mod(n, nPeriod)+1 > nBurst
%                         targets(n) = 0;
%                     end
%                 end
            else
                targets(nKnown+1:length(labels)) = 0;
            end
        end
    end
    
    methods        
        %labels should be of the form 0 = unknown, numbers > 0 refer to
        %class labels
        function obj = ExpRunner(expSet, synthSet, bccSet, labelledTrainingData,...
                labelledTestData, labels, basePost, agents, trainingLabels )
            obj.expSet = expSet;
            obj.bccSet = bccSet;
            obj.synthSet = synthSet;
            
            if exist('labels','var')
                obj.loadData = false;

                obj.labelledTrainingData = labelledTrainingData;
                obj.labelledTestData = labelledTestData;
                
                if size(labels, 1) > 1 && size(labels, 2)==1
                    labels = labels';
                end               
                
                obj.labels = labels;
                if ~isempty(labels)
                    obj.expSet.nSamples = length(obj.labels);                
                end
            end
            
            if exist('basePost','var')
                if iscell(basePost)
                    obj.basePost = [];
                    obj.baseLogodds = [];
                    obj.baseScoreSet = basePost;
                    
                    if isempty(labels)
                        obj.expSet.nSamples = max(obj.baseScoreSet{2});
                    end
                else
                    obj.basePost = basePost;
                    obj.baseLogodds = -log((1 ./ basePost) - 1);
                    obj.baseScoreSet = [];
                    
                    if isempty(labels)
                        obj.expSet.nSamples = length(basePost);
                    end
                end
            end
            
            if exist('trainingLabels','var')
                obj.trainingLabels = trainingLabels;
            end
            
            if exist('agents','var')
                obj.agents = agents;
            end
            
            if isempty(obj.labels)
                %obj.labels = -1 * ones(obj.expSet.nSamples, 1);
                obj.labels = sparse(double(obj.expSet.nSamples), 1);
            end
        end

        function [targets, nKnown] = getTrainingLabels(obj)
            if isempty(obj.trainingLabels)
                %No. target labels we know: usually the first nKnown labels.
                nKnown = round(obj.expSet.propKnown(obj.expSet.iPropKnown) ...
                            * obj.expSet.nSamples);
                if size(obj.labels, 1) > size(obj.labels, 2)
                    obj.labels = obj.labels';
                end
                if isempty(obj.labels)
                    display('no target labels');
                    targets = -1 * ones(1, obj.expSet.nSamples);
                else
                    targets = ExpRunner.getKnownTargets(...
                        obj.labels(1:obj.expSet.nSamples), nKnown, ...
                        obj.expSet.spreadTargets);
                end
                obj.trainingLabels = targets;
            else
                targets = obj.trainingLabels';
                nKnown = sum(targets~=0);
            end
        end
        
        function [nCombiners, nonDetIdx, targets] = createCombiners(obj, combMethods, nClusters, runNonDetOnly)
            
            [obj.combiner nonDetIdx targets] = CombinerFactory.createCombiners(...
                obj, combMethods, nClusters, runNonDetOnly, obj.bccSet.minScore, obj.bccSet.maxScore);
            nCombiners = length(obj.combiner);            
        end
        
        function [combinedPost, agentRatings, labels, basePost, combinedPostKnowns, labelsKnowns, baseKnowns] = runCombineAll(obj, ...
                newData, runBaseAgents, drawGraphs, combMethods, sortResults, runNonDetOnly)
            
            %train and test the agents (but no clustering)
            obj.runNoClustering(newData, runBaseAgents);  
            
            if nargin < 6
                sortResults = false;
            end
            
            if nargin < 7
                runNonDetOnly = false;
            end
            [combinedPost, agentRatings, combinedPostKnowns, labelsKnowns] = obj.combineAll(drawGraphs,...
                combMethods, sortResults, runNonDetOnly);
            
            basePost = obj.basePost;
            if ~isempty(basePost)
                baseKnowns = basePost(:, obj.getTrainingLabels()==0);
            end
            
            labels = obj.labels;
            
%             display('result from exprunner:');
%             for c=1:numel(combMethods)
%                 display(num2str(round(100*abs(labels(1:size(combinedPost,2))-combinedPost(c,:)-1))));
%             end
        end
                    
        function [combinedPost, agentRatings, combinedPostKnowns, labelsKnowns] = combineAll(obj, ...
                drawGraphs, combMethods, sortResults, runNonDetOnly)
            
            if ~exist('drawGraphs', 'var')
                drawGraphs = true;
            end
            if nargin < 5
                runNonDetOnly = false;
            end

            %find the best-performing agent for each data point. Use this
            %for comparison with realistic combination and selection
            %methods.
            
            [nCombiners, nonDetIdx, targets] = obj.createCombiners(combMethods, 1, runNonDetOnly);                
            
            if obj.expSet.nClasses==2
                combinedPost = zeros(nCombiners, obj.expSet.nSamples);
            else
                combinedPost = cell(1,obj.expSet.nClasses);
                for j=1:obj.expSet.nClasses
                    combinedPost{j} = zeros(nCombiners, obj.expSet.nSamples);
                end
            end
            agentRatings = cell(nCombiners, 1);     
            
            legendStrings = {nCombiners};
            
            combinerInfo = {nCombiners};
            
            if ~isempty(obj.baseScoreSet)
                combinerInputData = obj.baseScoreSet;
                if obj.convertTargetsToAgents
                    combinerInputData{4} = targets;
                end     
            else
                combinerInputData = obj.basePost(:, 1:obj.expSet.nSamples);        
                if obj.convertTargetsToAgents
                    combinerInputData = [combinerInputData; targets];
                end
            end
            
            if isempty(obj.runTime)
                obj.runTime = zeros(1, nCombiners);
            end
            
            %run combination functions over the agents' outputs.
            for i=1:nCombiners
                
                if runNonDetOnly && sum(nonDetIdx==i)==0
                    combinedPost(i, :) = -1;
                    continue;
                end
                
                legendStrings{i} = obj.combiner{i}.label;
                                 
                if ~iscell(combinerInputData) && size(combinerInputData, 1)<1000
                    obj.combiner{i}.useSparse = false;
                else
                    obj.combiner{i}.useSparse = true;
                end
                resultsPost = obj.combiner{i}.correctBaseData(combinerInputData);
                %resultsPost = combinerInputData;
                
                s=cputime;
                [combinedDecisions, agentRatings{i}] = ...
                    obj.combiner{i}.combineDecisions(resultsPost);
                t=cputime-s;
                
                obj.runTime(i) = obj.runTime(i) + t;
                if obj.expSet.nClasses ==2
                    combinedPost(i, 1:size(combinedDecisions,2)) = combinedDecisions;
                else
                    for j=1:obj.expSet.nClasses
                        if size(combinedDecisions, 1) == 1
                            decisions = combinedDecisions(1,:);
                        else
                            decisions = combinedDecisions(j,:);
                        end
                        combinedPost{j}(i, 1:size(combinedDecisions,2)) = decisions;
                    end
                end
                
                %normalise the scores
%                 if ~obj.combiner{i}.normalised
%                     combinedPost(i, :) = (combinedPost(i, :) - obj.bccSet.minScore) ./ ...
%                         (obj.bccSet.maxScore-obj.bccSet.minScore);
%                 end
                
                combinerInfo{i} = obj.combiner{i}.combinerInfo;
            end
            
            obj.combinedPost = combinedPost;
            
            if obj.saveData
                save(sprintf('%s%i_%i__%i_combinerInfo.mat', ...
                    obj.expSet.getDataDir(), ...
                    obj.expSet.iLambdaMag, ...
                    obj.expSet.iLambdaSym, ...
                    obj.expSet.iNu), 'combinerInfo');
            end
            %sort results in order of the first combination - to make
            %nicer-looking graphs
            if sortResults
                if obj.expSet.nClasses<=2
                    combinedPost = sortrows(combinedPost')';
                else
                    for j=1:obj.expSet.nClasses
                        combinedPost{j} = sortrows(combinerPost{j}')';
                    end
                end
            end

%             if we have no true labels, use the first combination function
            %for error/comparison
            if isempty(obj.labels)
                obj.labels = combinedPost(1, :);
            end
            
            if obj.expSet.nClasses<=2
                combinedPostKnowns = combinedPost(:, obj.getTrainingLabels()==0);
            else
                combinedPostKnowns = {obj.expSet.nClasses};
                for j=1:obj.expSet.nClasses
                    combinedPostKnowns{j} = combinedPost{j}(:, obj.getTrainingLabels()==0);
                end
            end
            labelsKnowns = obj.labels(obj.getTrainingLabels()==0);
            
            if drawGraphs
                legendStrings{length(legendStrings)+1} = 'real answers';
                
                %plot only data for which we have the correct labels and
                %where labels were not given in training
                if isempty(obj.changes) && isempty(obj.expSet.definiteSwitch)
                    testLabels = obj.labels~=0 & obj.getTrainingLabels()==0;
                else
                    if isempty(obj.changes)
                        changeAt = round(obj.expSet.definiteSwitch.*obj.expSet.nSamples);
                        obj.changes = cell(1,4);
                        obj.changes{1} = changeAt;
                        obj.changes{2} = ones(1,length(changeAt));
                        obj.changes{3} = obj.changes{1};
                        obj.changes{4} = obj.changes{2};
                    end
                    testLabels = obj.labels~=0;
                    display('As we are showing are changepoints, and the training data is mixed with the test data, we show the training data responses as well as test data responses.');
                end
                
                if obj.expSet.nClasses==2
                    labelledCombinedPost = combinedPost(:, testLabels);
                    knownLabels = obj.labels(testLabels);

                    lcPostWithLabels = [labelledCombinedPost; knownLabels-1];

                    if sortResults
                        [sortedMeans, idxs] = sort(lcPostWithLabels(1,:));
                        lcPostWithLabels = lcPostWithLabels(:, idxs);
                    end
                    graph = graphs.ClassifierPerformanceGraph(obj.expSet, ...
                        {}, 1, legendStrings, obj.changes);

    %                 graph.drawResultHistogram(lcPostWithLabels, 'Histogram of scores', {'not SN' 'SN'});
                    graph.drawPostGraph(lcPostWithLabels, 'Output of Combination Methods');
                    graph.drawErrorGraph(knownLabels, labelledCombinedPost, 'Absolute Error of Combination Methods');
                else
                    for j=1:obj.expSet.nClasses-1
                        labelledCombinedPost = combinedPost{j}(:, testLabels);
                        knownLabels = obj.labels(testLabels);

                        lcPostWithLabels = [labelledCombinedPost; knownLabels-1];

                        if sortResults
                            [sortedMeans, idxs] = sort(lcPostWithLabels(1,:));
                            lcPostWithLabels = lcPostWithLabels(:, idxs);
                        end
                        graph = graphs.ClassifierPerformanceGraph(obj.expSet, ...
                            {}, 1, legendStrings, obj.changes);

        %                 graph.drawResultHistogram(lcPostWithLabels, 'Histogram of scores', {'not SN' 'SN'});
                        graph.drawPostGraph(lcPostWithLabels, 'Output of Combination Methods');
                        graph.drawErrorGraph(knownLabels, labelledCombinedPost, 'Absolute Error of Combination Methods');

                    end
                end
%                 knownLabels = knownLabels-1;
%                 graphs.ClassifierPerformanceGraph.drawRoc(labelledCombinedPost, knownLabels, combMethods);
            end
            if obj.saveData
                save(sprintf('%scombineAll.mat', obj.expSet.getDataDir()));
            end
        end

        function runCombineByCluster(obj, newData, runBaseAgents, combMethods, methodInDetail)
            
            %train and test the agents (but no clustering)
            obj.runWithClustering(newData, runBaseAgents);  
            
            K = obj.clusterer.K;
            nCombiners = obj.createCombiners(combMethods, K);
            
            nResults = nCombiners*K + nCombiners; %add an extra set for when we ignore the clustering
            
            obj.combinedPost = zeros(nResults, obj.expSet.nSamples);
            combinedLogodds = zeros(nResults, obj.expSet.nSamples);            
            
            K = obj.clusterer.K;
            agentClusterMembership = obj.clusterer.agentClusterMembership;
            
            legendStrings = {nCombiners*K};
            
            for i=1:nCombiners
                rowIdxs = (1:K) .* nCombiners - nCombiners + i;
                
                %run combination functions over the agents' outputs.
                [obj.combinedPost(rowIdxs, :), combinedLogodds(rowIdxs, :)] = ...
                    obj.combiner{i}.combineDecisions(obj.basePost, ...
                        obj.baseLogodds, agentClusterMembership);
            end
            
            %get the combinations using no clustering.
            [obj.combinedPost(nCombiners*K+1:nCombiners*K+nCombiners, :)] =...
                obj.combineAll(false, combMethods);
            
            %draw the graphs with a window for each cluster, so we can
            %compare each combination method for a particular cluster.
            %May want to alter graph_agent_performance so we can put the
            %graphs for all methods into a single subplot - not sure if
            %this will improve or hinder legibility.
            
            if exist('methodInDetail','var')
                for k=1:K
                    legendStrings{k} = sprintf('%s, cluster %d', obj.combiner{2}.label, k);
                end

                legendStrings{K+1} = sprintf('%s, all agents', obj.combiner{2}.label);
            
                compareIdx = methodInDetail + nCombiners*K;
                %show the graphs that compare the XXX combination method for
                %each cluster to the second combination method for all agents.
                for k=1:K
                    graph = graphs.ClassifierPerformanceGraph(obj.expSet, ...
                        {}, 2, {legendStrings{k} legendStrings{K+1}}, obj.changes); 

                    graph.drawErrorGraph(obj.labels, obj.combinedPost([methodInDetail+nCombiners*(k-1) compareIdx], :), ...
                        sprintf('%s for cluster %i', obj.combiner{methodInDetail}.label, k));
                end
            end
            
            for k=1:K+1
                
                for i=1:nCombiners
                    legendStrings{i} = sprintf('%s, cluster %d', obj.combiner{i}.label, k);
                end

                graph = graphs.ClassifierPerformanceGraph(obj.expSet, ...
                    {}, 1, legendStrings, obj.changes);
            
                firstRow = k*nCombiners - nCombiners + 1;
                lastRow = k*nCombiners;
                
                graph.drawErrorGraph(obj.labels, obj.combinedPost(firstRow:lastRow, :), ...
                    sprintf('combination methods for cluster %i', k));
                
                graph.drawPostGraph(obj.combinedPost(firstRow:lastRow, :), ... 
                    sprintf('combination methods - posterior probability for cluster %i', k));
            end
            graph = graphs.ClassifierPerformanceGraph(obj.expSet, ...
                    {}, 1, legendStrings, obj.changes, obj.agents);
            graph.drawSensorSubscriptions();
            if obj.saveData
                save(sprintf('%scombine_by_cluster.mat', obj.expSet.getDataDir()));
            end
        end        
        
        %functions related to running base classifiers
        
        function loadDataForAgents(obj, testDataFile, trainingDataFile)
            obj.labelledTestData = dlmread(testDataFile);
            obj.labelledTrainingData = dlmread(trainingDataFile);
            obj.labels = obj.labelledTestData( ...
                :, ...
                obj.expSet.nSensors()+1)';

            obj.changes = datageneration.sensorChangepoints(...
                                obj.expSet, obj.labelledTestData);
        end
        
        function generateDataForAgents(obj)
            [obj.labelledTrainingData obj.labelledTestData obj.labels] = ...
                datageneration.generateData(obj.synthSet, obj.expSet.expLabel, obj.expSet.nSamples);

            obj.changes = ...
                datageneration.sensorChangepoints(...
                                obj.synthSet, obj.expSet.nSamples, obj.labelledTestData);           
        end
        
        function [errors fnr_agents fpr_agents] = agentErrorRates(obj)
            alabels = repmat(obj.labels-1, length(obj.agents), 1);
            errors = sum( (obj.basePost-alabels).^2, 2) ./ length(obj.labels);
            
            roundedDecs = round(obj.basePost);
            posIdxs = obj.labels==2;
            negIdxs = obj.labels==1;
            fnr_agents = sum(1-obj.basePost(:, posIdxs),2) ./ sum(posIdxs);
            fpr_agents = sum(obj.basePost(:,negIdxs),2) ./ sum(negIdxs);
        end
        
        function runAgents(obj, agentsFile)
            %ensure the base agents code can be found
            if strcmp(obj.synthSet.agentType, 'dlr')
                addpath(obj.synthSet.dlrPath);
            end

            [obj.agents obj.basePost obj.baseLogodds obj.agentConfidence] = ...
                agentexecution.trainAndTestRun( ...
                    obj.synthSet, ...
                    obj.labelledTrainingData, ...
                    obj.labelledTestData, ...
                    obj.synthSet.reuseSubs ...
                );

            if obj.synthSet.flipAgents > 0
                if obj.synthSet.nInformedAgents > 0
                    nFlips = obj.synthSet.flipAgents * obj.synthSet.nInformedAgents;
                else
                    nFlips = obj.synthSet.flipAgents * obj.synthSet.nAgents;
                end
                
                for a=1:nFlips
                    obj.basePost(a, :) = 1 - obj.basePost(a, :);
                    obj.baseLogodds(a, :) = - obj.baseLogodds(a, :);
                end
            end

            save(agentsFile, 'obj');             
        end
        
        function runFakeAgents(obj, agentsFile)
            %agent map should contain a 3D matrix where layers correspond to
            %agents, rows correspond to the correct label class,
            % and column entries are probability of getting each class
            % label given the correct class
            agentMap = obj.synthSet.fakeAgents;
            
            results = zeros(obj.synthSet.nAgents, obj.expSet.nSamples);
            for n=1:length(obj.labels)
                for a=1:obj.synthSet.nAgents
                    cdf = cumsum( agentMap(obj.labels(n), :, a) );
                    results(a, n) = sum(cdf < rand*cdf(end)); %labels start at 0
                end
            end
            
            if ~isempty(obj.synthSet.definiteSwitch)
                
                obj.changes = cell(1,4);
                
                for sw=round(obj.synthSet.definiteSwitch*obj.expSet.nSamples)
                    order = randperm(obj.synthSet.nAgents);
                    display(['shuffling agents: ' num2str(order)]);
                    while sum(order==(1:obj.synthSet.nAgents))==obj.synthSet.nAgents
                        order = randperm(obj.synthSet.nAgents);
                    end
                    oldResults = results;
                    for a=1:obj.synthSet.nAgents
                        results(order(a), sw:end) = oldResults(a, sw:end);
                    
                    
                        obj.changes{1} = [obj.changes{1} sw];
                        obj.changes{3} = obj.changes{1};

                        obj.changes{2} = [obj.changes{2} order(a)];
                        obj.changes{4} = [obj.changes{4} a];
                    end
                end                
            end
            
            obj.basePost = results;
            
            datageneration.checkDataDir(obj.expSet);
            save(agentsFile, 'obj');   
        end
        
        function loadAgentsResults(obj, agentsFile)
            savedRunner = load(agentsFile, 'obj');

            obj.basePost = savedRunner.obj.basePost;
            obj.baseLogodds = savedRunner.obj.baseLogodds;
            obj.agents = savedRunner.obj.agents;
            obj.agentConfidence = savedRunner.obj.agentConfidence;
        end
        
        function runNoClustering(obj, newData, runBaseAgents)
            
            if strcmp(runBaseAgents, 'dontTouch')
                display('Not running or reloading agents.');
                return;
            end
            
            %data for base classifiers ------------------------------------
            testDataFile = sprintf('%s%s_test_data.mat', ...
                obj.expSet.getDataDir(), ...
                obj.expSet.expLabel);
            
            trainingDataFile = sprintf('%s%s_training_data.mat', ...
                obj.expSet.getDataDir(), ...
                obj.expSet.expLabel);    
            
            if ~exist(testDataFile, 'file') && ~exist(trainingDataFile, 'file')
                if isempty(find(obj.labels, 1)) && ~strcmp(runBaseAgents, 'loadSaved')
                    newData = true;
                end
            end
            
            labelsFile = sprintf('%s%s_labels.mat', ...
                obj.expSet.getDataDir(), ...
                obj.expSet.expLabel);
            
            if exist(labelsFile, 'file') && ~strcmp(runBaseAgents, 'newAgents') ...
                    && (isempty(obj.labels) || strcmp(runBaseAgents, 'loadSaved'))
%                newData = false;
                obj.loadData = false;
                obj.labels = dlmread(labelsFile);
            end
            
            if newData
                if strcmp(runBaseAgents, 'loadSaved')
                    display('You tried to create new data or the old labels were missing')
                    display('-- but you will not create new agents. Quite frankly sir, you are mad.');
                elseif strcmp(runBaseAgents, 'newFakeAgents')
                    display('Creating labels for fake agents');
                    obj.labels = datageneration.generateLabels(2, obj.expSet.nSamples, obj.expSet);
                else
                    obj.generateDataForAgents();                    
                end                 
            elseif obj.loadData
                obj.loadDataForAgents(testDataFile, trainingDataFile);
            else
                display('Not touching data for simulated base classifiers');
            end
            
            if max(obj.labels) < obj.expSet.nClasses && min(obj.labels)==0
                obj.labels = obj.labels + 1;
            end

            agentsFile = sprintf('%sagents.mat', ...
                    obj.expSet.getDataDir());   
                
            % load, rerun or ignore base classifiers/agents----------------
            if strcmp(runBaseAgents, 'newAgents') || ...
                    (strcmp(runBaseAgents, 'loadSaved') && ~exist(agentsFile, 'file'))
                obj.runAgents(agentsFile);
            elseif strcmp(runBaseAgents, 'newFakeAgents')
                display('creating fake agents and their results!');
                obj.runFakeAgents(agentsFile);
            elseif strcmp(runBaseAgents, 'loadSaved')
                obj.loadAgentsResults(agentsFile);
            end
        end
        
        function runWithClustering(obj, newData, runBaseAgents)

            obj.runNoClustering(newData, runBaseAgents);

            obj.clusterer = clustering.agentClustererKmeans(obj.expSet);
            obj.clusterer.cluster(obj.basePost, obj.baseLogodds, ...
                obj.labels, obj.agentConfidence);
            if obj.saveData
                save(sprintf('%swith_clustering.mat', obj.expSet.getDataDir()));
            end
        end
            
        function runWithGraphs(obj, newData, runBaseAgents, filePrefix, graphVisible)

            if exist('filePrefix', 'var')
                obj.expSet.expLabel = filePrefix;        
            end
            
            if exist('graphVisible', 'var')
                obj.expSet.graphVisible = graphVisible;
            end        
        
            obj.runWithClustering(newData, runBaseAgents);
        
            perfGraph = graphs.ClassifierPerformanceGraph(obj.expSet, ...
                obj.clusterer.agentClusterMembership, 0, {}, obj.changes, obj.agents);
            
            perfGraph.drawGraphs(obj.basePost, obj.baseLogodds, obj.labels, 'individualAgents');
            
            clusterGraph = graphs.ClusterGraph(obj.expSet, obj.changes, obj.clusterer);
            
            clusterGraph.drawGraphs(obj.basePost);
            if obj.saveData
                save(sprintf('%swith_graphs.mat', obj.expSet.getDataDir()));
            end
        end
    end
end