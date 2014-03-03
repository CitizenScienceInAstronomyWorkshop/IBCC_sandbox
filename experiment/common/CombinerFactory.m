classdef CombinerFactory

    methods (Static)

        function [combiner] = createCombiner(runner, combMethod, targets, ...
                bccTargets, nClusters, nAgents, nKnown, minScore, maxScore, windowSize, seq)
            
            %some methods learn froma window of data only
            if ~exist('windowSize','var')
                windowSize = length(targets);
            end
            
            %run sequential methods in batch or sequential mode?
            if ~exist('seq','var')
                seq = true;
            end
            
            expSet = runner.expSet;
            bccSet = runner.bccSet;            
            
            if strcmp(combMethod, combiners.TrueBestAgent.shortLabel)
                display('Known labels used incorrectly for TrueBestAgent');
                combiner = combiners.TrueBestAgent(...
                    runner.labels, nAgents, nClusters);

            %Avoid weight calculations; reject poorer classifiers; used
            %with limited resources or when only one base classifier is
            %informed.
            elseif strcmp(combMethod, combiners.SelectBest.shortLabel)
                combiner = combiners.SelectBest(...
                    targets, nAgents, 20, 10, nClusters, bccSet.targetsAsAgents);

            %Simple majority voting; use voting simplifies calculations
            %Doesn't give confidence level; reduces effect of noise
            elseif strcmp(combMethod, combiners.SimpleMajorityVoting.shortLabel)
                mv = combiners.SimpleMajorityVoting(nAgents, nClusters,expSet.nClasses);
                mv.noScore = expSet.noScore;
                mv.minScore = minScore;
                mv.maxScore = maxScore;
                mv.voteThreshold = expSet.voteThreshold;
                combiner = mv;
            elseif strcmp(combMethod, combiners.SimpleMajorityVoting.shortLabelUnrounded)
                mv = combiners.SimpleMajorityVoting(nAgents, nClusters);
                mv.roundCombination = false;
                mv.noScore = expSet.noScore;
                mv.minScore = minScore;
                mv.maxScore = maxScore;
                mv.voteThreshold = expSet.voteThreshold;
                combiner = mv;

            %Simple mean; assume uncorrelated errors, similar confidence/prob & size
            %of error in each base learner
            elseif strcmp(combMethod, combiners.MeanDecision.shortLabel)
                mean = combiners.MeanDecision(nAgents, nClusters, ...
                    expSet.nClasses, maxScore, minScore);
                mean.noScore = expSet.noScore;
                combiner = mean;
            elseif strcmp(combMethod, combiners.RoundedMeanDecision.subShortLabel)
                mean = combiners.RoundedMeanDecision(nAgents, nClusters);
                mean.noScore = expSet.noScore;
                combiner = mean;
            elseif strcmp(combMethod, combiners.CountDecisions.shortLabel)
                combiner = combiners.CountDecisions(nAgents, nClusters);
                combiner.noScore = expSet.noScore;
                
            %DLR; treat base classifier outputs as features
            elseif strcmp(combMethod, combiners.DlrCombiner.shortLabel)
                dlrComb = combiners.DlrCombiner(nAgents, nClusters, ...
                    targets, runner.expSet.dlrPath);
                dlrComb.noScore = expSet.noScore;
                combiner = dlrComb;
            elseif strcmp(combMethod, combiners.DlrCombiner.shortLabelRounded)
                dlrComb = combiners.DlrCombiner(nAgents, nClusters, ...
                    targets, runner.expSet.dlrPath);  
                dlrComb.roundCombination = true;
                dlrComb.noScore = expSet.noScore;
                combiner = dlrComb;
            elseif strcmp(combMethod, combiners.LrCombiner.shortLabel)
                lrComb = combiners.LrCombiner(nAgents, nClusters, targets);
                combiner = lrComb;
                
            %NB estimates of logodds: recalculates the precision of 
            %each base learner in order to combine their estimates     
            elseif strcmp(combMethod, combiners.weighted.NaiveBayesLogodds.shortLabel)
                weighted = combiners.weighted.NaiveBayesLogodds(nAgents, ...
                    nClusters, targets, runner.agents, 20, ...
                    runner.labelledTestData(:, 1:runner.synthSet.nSensors()));
                combiner = weighted;
            %NB estimates of logodds; 
            elseif strcmp(combMethod, combiners.weighted.NaiveBayesLogodds.precalcShortLabel)
                if isempty(runner.labelledTestData)
                    weighted = combiners.weighted.NaiveBayesLogodds(nAgents, ...
                    nClusters, targets, runner.agents, 20); 
                    combiner = weighted;
                else
                    data = runner.labelledTestData(:, 1:runner.synthSet.nSensors());
                    weighted = combiners.weighted.NaiveBayesLogodds(nAgents, ...
                        nClusters, targets, runner.agents, 20, data); 
                    combiner = weighted;
                end
            elseif strcmp(combMethod, combiners.weighted.NaiveBayesLogodds.precalcShortLabelRounded)
                weighted = combiners.weighted.NaiveBayesLogodds(nAgents, ...
                    nClusters, targets, runner.agents, 20); 
                weighted.roundCombination = true;
                combiner = weighted;

            %Weighted Sum: when the combination makes a mistake, 
            %weights of base learners are adjusted proportional to
            %their contribution to the incorrect decision
            elseif strcmp(combMethod, combiners.weighted.WeightedSum.shortLabel)
                weighted = combiners.weighted.WeightedSum(nAgents, ...
                    nClusters, targets, runner.agents);
                weighted.onTheFly = false;
                weighted.noScore = expSet.noScore;
            elseif strcmp(combMethod, combiners.weighted.WeightedSum.shortLabelRounded)
                weighted = combiners.weighted.WeightedSum(nAgents, ...
                    nClusters, targets, runner.agents);
                weighted.onTheFly = false;
                weighted.roundCombination = true;
                weighted.noScore = expSet.noScore;

            %Weighted Majority: similar weights to weighted sum
            elseif strcmp(combMethod, combiners.weighted.WeightedMajority.subShortLabel)
                weighted = combiners.weighted.WeightedMajority(nAgents, ...
                    nClusters, targets, runner.agents);
                weighted.onTheFly = false;
                weighted.postScores = expSet.postScores;
                weighted.voteThreshold = expSet.voteThreshold;
                weighted.nVotesPerClassifier = 3;
                weighted.noScore = expSet.noScore;
            elseif strcmp(combMethod, combiners.weighted.WeightedMajority.subShortLabelUnrounded)
                weighted = combiners.weighted.WeightedMajority(nAgents, ...
                    nClusters, targets, runner.agents);
                weighted.onTheFly = false;
                weighted.roundCombination = false;
                weighted.voteThreshold = expSet.voteThreshold;                    

            %Bma
            elseif strcmp(combMethod, combiners.weighted.Bma.subSubShortLabel)
                weighted = combiners.weighted.Bma(nAgents, nClusters, targets, ...
                    runner.agents, windowSize);
                weighted.onTheFly = false;
                combiner = weighted;

            %Ibcc
            %Likely to perform well on static problems only at the
            %moment, as confusion matrix doesn't account for changes
            %over time
            elseif strcmp(combMethod, combiners.bcc.IbccSampling.shortLabel)
                ibcc = combiners.bcc.IbccSampling(bccSet, nAgents, nClusters, bccTargets, ...
                    runner.agents, expSet.nClasses, expSet.nScores);
                %use a sampled Alpha or Alpha0+voteCounts
                ibcc.fixedAlpha = true;
            elseif strcmp(combMethod, combiners.bcc.IbccVb.shortLabelSeq)
                ibcc = combiners.bcc.IbccVb(bccSet, nAgents, nClusters, bccTargets, ...
                    runner.agents, expSet.nClasses, expSet.nScores);
                ibcc.sequential = true;
                ibcc.label = [ibcc.label ' sequential'];               
            elseif strcmp(combMethod, combiners.bcc.IbccVb.shortLabel)
                ibcc = combiners.bcc.IbccVb(bccSet, nAgents, nClusters, bccTargets, ...
                    runner.agents, expSet.nClasses, expSet.nScores);
                ibcc.sequential = false;                   
            elseif strcmp(combMethod, combiners.bcc.IbccPooledVb.shortLabel)
                ibcc = combiners.bcc.IbccPooledVb(bccSet, nAgents, nClusters, bccTargets, ...
                    runner.agents, expSet.nClasses, expSet.nScores);
                ibcc.sequential = false;    
            elseif strcmp(combMethod, combiners.bcc.IbccDiff.shortLabel)
                ibcc = combiners.bcc.IbccDiff(bccSet, nAgents, nClusters, bccTargets, ...
                    runner.agents, expSet.nClasses, expSet.nScores);
                ibcc.sequential = false;                    
	    elseif strcmp(combMethod, combiners.bcc.DynIbccVb.shortLabel)
                ibcc = combiners.bcc.DynIbccVb(bccSet, nAgents, nClusters, bccTargets, ...
                    runner.agents, expSet.nClasses, expSet.nScores);
                ibcc.sequential = false;        

            elseif strcmp(combMethod, combiners.bcc.DynIbccVb.shortLabelSep)
                ibcc = combiners.bcc.DynIbccVb(bccSet, nAgents, nClusters, bccTargets, ...
                    runner.agents, expSet.nClasses, expSet.nScores);
                ibcc.sequential = false;        
                ibcc.sepConfMatrix = true;

            elseif strcmp(combMethod, combiners.bcc.DynIbccVb.shortLabelInitAll)
                ibcc = combiners.bcc.DynIbccVb(bccSet, nAgents, nClusters, bccTargets, ...
                    runner.agents, expSet.nClasses, expSet.nScores);
                ibcc.sequential = false;        
                ibcc.initAll = true;

            elseif strcmp(combMethod, combiners.bcc.ibccVbInit.IbccVbInitMeds.shortLabel)
                ibcc = combiners.bcc.IbccVb(bccSet, nAgents, nClusters, bccTargets, ...
                    runner.agents, expSet.nClasses, expSet.nScores);
                ibcc.initMeds = true;                
            elseif strcmp(combMethod, combiners.bcc.ibccVbInit.IbccVbInitAllLabs.shortLabel)
                ibcc = combiners.bcc.IbccVb(bccSet, nAgents, nClusters, bccTargets, ...
                    runner.agents, expSet.nClasses, expSet.nScores);
                ibcc.initUseAllLabs = true;             
            elseif strcmp(combMethod, combiners.bcc.CIbccVb.shortLabel)
                ibcc = combiners.bcc.CIbccVb(bccSet, nAgents, nClusters, bccTargets, ...
                    runner.agents, expSet.nClasses, expSet.nScores);                    
            elseif strcmp(combMethod, combiners.bcc.IbccVbAux.shortLabel)
                ibcc = combiners.bcc.IbccVbAux(bccSet, nAgents, nClusters, bccTargets, ...
                    runner.agents, expSet.nClasses, expSet.nScores);                     
            elseif strcmp(combMethod, combiners.bcc.IbccVbEmInit.shortLabel)
                ibcc = combiners.bcc.IbccVbEmInit(bccSet, nAgents, nClusters, bccTargets, ...
                    runner.agents, expSet.nClasses, expSet.nScores);
            elseif strcmp(combMethod, combiners.bcc.IbccEm.shortLabel)
                ibcc = combiners.bcc.IbccEm(bccSet, nAgents, nClusters, targets, ...
                    runner.agents, expSet.nClasses, expSet.nScores);
                display('WARNING: ICC EM may try to incorrectly use trustFinalAgent');
            elseif strcmp(combMethod, combiners.bcc.IbpcVb.shortLabel)
                ibcc = combiners.bcc.IbpcVb(bccSet, nAgents, nClusters, bccTargets, ...
                    runner.agents, expSet.nClasses, expSet.nScores,...
                    [bccSet.lambdaMag(bccSet.iLambdaMag) bccSet.trustedLambdaMag], ...
                    [bccSet.lambdaSym(bccSet.iLambdaSym) bccSet.trustedLambdaSym]);
            elseif strcmp(combMethod, combiners.bcc.IbpcVb2.shortLabel)
                ibcc = combiners.bcc.IbpcVb2(bccSet, nAgents, nClusters, bccTargets, ...
                    runner.agents, expSet.nClasses, expSet.nScores,...
                     [bccSet.lambdaMag(bccSet.iLambdaMag) bccSet.trustedLambdaMag], ...
                    [bccSet.lambdaSym(bccSet.iLambdaSym) bccSet.trustedLambdaSym]); 
            elseif strcmp(combMethod, combiners.bcc.CIbcc.shortLabel)
                combiner = combiners.bcc.CIbcc(bccSet, nAgents, nClusters, ...
                    runner.labels, runner.agents);

            %methods we are not testing ------------------------------

            %Like NB estimates of logodds, but uses bayes error to
            %calculate weights rather than true precision - should
            %remove?
            elseif strcmp(combMethod, combiners.weighted.ExpBayesErrorLogodds.subShortLabel)
                combiner = combiners.weighted.ExpBayesErrorLogodds(...
                    nAgents, nClusters, targets, runner.agents, 20);

            %Precision-Weighted Average/Sum; poorly motivated - should remove?
            elseif strcmp(combMethod, combiners.weighted.PrecWeightedAveraging.subShortLabel)
                combiner = combiners.weighted.PrecWeightedAveraging(...
                    nAgents, nClusters, runner.labels, runner.agents, 20);
            end
            
            if exist('ibcc', 'var')
                %extra settings for IBCC
                ibcc.nKnown = nKnown;
%                 ibcc.confSave = expSet.getDataDir();
%                 ibcc.datasetLabel = runner.datasetLabel;
                ibcc.noScore = expSet.noScore;
                                
                if ~isempty(bccSet.trustedAlpha)
                    ibcc.trustedAlpha = bccSet.trustedAlpha;
                end      
                ibcc.setTrustedAgents(bccSet.targetsAsAgents);   
                                
                if ~isempty(bccSet.IbccMapFunction)
                    ibcc.mapFunction = bccSet.IbccMapFunction;
                end
                
                combiner = ibcc;
            end 
            
            if exist('weighted', 'var')
                weighted.minScore = minScore;
                weighted.maxScore = maxScore;
                weighted.setSequential(seq);
                weighted.windowSize = windowSize;
%                 weighted.useFinalWeights = false;
                combiner = weighted;
            end
            
            combiner.setProgressMonitor(runner.progressMonitor);              
        end
        
        function [combiner nonDetIdx targets] = createCombiners(runner, ...
                combMethods, nClusters, runNonDetOnly, minScore, maxScore)
        %CombinerFactory Factory to create a set of combiners given a set of combination method strings 
        %string and the settings encapsulated in runner, which should be of type
        %ExpRunner.

            if isempty(runner.agents)
                if ~isempty(runner.baseScoreSet)
                    nAgents = max(runner.baseScoreSet{1});
                else
                    nAgents = size(runner.basePost, 1);
                end
            else
                nAgents = length(runner.agents);            
            end

            if runner.convertTargetsToAgents
                nAgents = nAgents + length(runner.targetsAsAgents);
            end

            nCombiners = length(combMethods);
            combiner = cell(nCombiners);

            if nargin < 4
                runNonDetOnly = false;
            end
            nonDetIdx = [];

            [targets, nKnown] = runner.getTrainingLabels();
            
%             if runner.targetsAsAgents
%                 %get no known targets for BCC - targets will be passed in
%                 %differently
%                 bccTargets = ExpRunner.getKnownTargets(...
%                     runner.labels(1:bccSet.nSamples), 0, ...
%                     bccSet.spreadTargets);
%             else
                bccTargets = targets;
%             end

            for i=1:length(combMethods)
                [combiner{i}] = CombinerFactory.createCombiner...
                    (runner, combMethods{i}, targets, bccTargets, ...
                    nClusters, nAgents, nKnown, minScore, maxScore);
                
                combiner{i}.setId(i);
                
                if runNonDetOnly && combiner{i}.nonDet
                    nonDetIdx = [nonDetIdx i];
                end 
            end
        end

    end
end