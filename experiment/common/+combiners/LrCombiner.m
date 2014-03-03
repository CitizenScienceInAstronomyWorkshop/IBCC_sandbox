classdef  LrCombiner < combiners.AbstractCombiner
    %Logistic regression combiner
    
    properties
        label = 'Logistic Regressor';
    end
    properties (Constant)
        shortLabel = 'lr';
    end
    
    methods
        function obj = LrCombiner(nAgents, K, targets)
            obj@combiners.AbstractCombiner(nAgents, K, targets);
        end
        
        function [combinedPost, agentRatings] = combineDecisions(obj, baseOutputs, ~)
               
            nClasses = length(unique(obj.targets)) - 1;
            if sum(obj.targets==-1)==0
                nClasses = nClasses + 1;
            end
               
            if nClasses > 2
                combinedPost = zeros(nClasses, length(baseOutputs));
            else
                combinedPost = zeros(1, length(baseOutputs));
            end
            
            T = obj.targets;
            if length(T) < size(baseOutputs, 2)
                T = [T -1*ones(1, size(baseOutputs, 2)-length(T))]';
            else
                T = T(1:size(baseOutputs,2));
            end

            Ttr = full(T(T~=-1)') + 1;
            Dtr = full(baseOutputs(:, T~=-1)') + 1;
            Dtest = full(baseOutputs(:,T==-1)') + 1; 
            
            %standard naive bayes
            B = NaiveBayes.fit(Dtr, Ttr, 'Distribution', 'mn');
            p = B.posterior(Dtest);
            
            combinedPost(:, T==-1) = p';
            for j=1:nClasses
                combinedPost(j, T~=-1) = Ttr;
            end
            agentRatings = B;
            
%             %naive bayes logop for each class -- too much memory; not
%             sure it makes sense anyway for discrete base classifiers?
%             for j=1:nClasses
%                 nbc = combiners.weighted.NaiveBayesLogodds(obj.nAgents, 1, obj.targets==j, [], size(baseOutputs,2), baseOutputs);
%                 combinedPost(j,:) = nbc.combineDecisions(baseOutputs);
%             end     
                        
            %             nBase = size(baseOutputs,1);
%             sensors = ones(1,nBase);
%             
%             nTargets = length(obj.targets);
%                         
%             lrs = {nClasses-1};
%             

%             %training
%             if nClasses==2
%                 start = 2;
%             else
%                 start = 1;
%             end
%             for j=start:nClasses
%                 
%                
%                 
%                 Ttr = sparse(Ttr==(j-1))';
%                 
%                 lrs{j} = baselearners.logisticRegressionAgent(sensors, 'lrc', nTargets, 1, 0);
%                 lrs{j}.train(Dtr, Ttr);
%                 
%                 Ttest = T(T==-1)';
%                                
%                 testResults = lrs{j}.test(Dtest, Ttest);
%                 combinedPost(j, T==-1) = testResults;
%                 combinedPost(j, T~=-1) = Ttr;
%             end
        end
    end
    
end

