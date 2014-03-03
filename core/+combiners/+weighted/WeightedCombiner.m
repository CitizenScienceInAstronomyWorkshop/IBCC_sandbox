classdef  WeightedCombiner < combiners.AbstractCombiner
    
    properties
        
        %prior for probability of class 1 
        %(in two class problem with classes 0 and 1)
        p_c1
        
        agents
        
        windowSize
                
        weightCalculator
        
        scoreSet = false
   
        sequential = false;
        
        useFinalWeights;
        
        maxScore = 1; %if these are different, we need to adjust so we can calculate (uncertain) votes.
        minScore = 0;        
    end
    
    methods (Abstract)
                
        combinedPost = combinationStep(obj, baseData, W, n);
        W = createWeightCalculator(obj, members, nSamples); 
    end
    
    methods
        function obj = WeightedCombiner(nAgents, K, targets, agents, windowSize)
            obj@combiners.AbstractCombiner(nAgents, K, targets);
            obj.agents = agents;
            if ~exist('windowSize', 'var')
                windowSize = 100000;
            end
            
            obj.windowSize = windowSize;
        end 
        
        function setSequential(obj, seq)
            obj.sequential = seq;
            obj.useFinalWeights = seq;
        end  
        
        function combinedPost = combineCluster(obj, baseData, T, members)
            
            %this may need to know which agents are in the cluster!
            if nargin < 4
                members = [];
            end
            obj.createWeightCalculator(members, size(T, 2));
           
            %if the weight calculation depends on the
            %combination step you will need to overwrite this - see
            %weighted sum. Can we organise these classes better?
            W = obj.weightCalculator.calculateWeights(baseData, baseData, T);

            combinedPost = obj.combinationStep(baseData, W);
        end
        
        function [combinedPost, W] = combineDecisions(obj, baseOutputs, clusters)
                        
            W = [];
            
            if obj.scoreSet
                nSamples = max(baseOutputs{2});
            else
                nSamples = size(baseOutputs, 2);
            end
            
            combinedPost = zeros(obj.K, nSamples);
            combinedLogodds = zeros(obj.K, nSamples);
                        
            if nargin < 3
                clusters = [];
            end
            
            %get the relevant inputs
            if obj.K == 1
                combinedPost(1, :) = ...
                    obj.combineCluster( baseOutputs, ...
                        obj.targets, clusters);    
            else
                for k=1:obj.K
                    
                    %clusters change over time, so dlr needs to take different
                    %inputs over time - this is not possible.
                    %So we set values for base learners outside the cluster
                    %to 0.5 (regardless of label) as we can't ignore them.
                   
                    inputPosts = baseOutputs;
                    inputPosts(clusters~=k) = 0.5;
                                        
                    [combinedPost(k, :), combinedLogodds(k, :)] = ...
                        obj.combineCluster( inputLogodds, inputPosts, ...
                            obj.targets, (clusters==k), k);
                end
            end
            
        end
    end
    
end

