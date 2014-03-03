classdef AbstractCombiner < matlab.mixin.Copyable
    %ABSTRACTCOMBINER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = protected)
        targets
        nAgents
        
        %number of clusters. We combine by cluster. If the clusters are not
        %to be used just set K = 1.
        K

        %progress monitor is an object with method monitorProgress(results)
        %that is called iteratively to measure the progress of combiners
        %that work through a number of steps. Used to measure how quickly
        %the combiner reaches a satisfactory solution.        
        progressMonitor = [];
        
        %number that can be assigned to an instance of a combiner for
        %reference
        id;
        
        nClasses;
        
        debug = true; %show debug/progress info or not
    end
    
    properties
        useSparse = true; %sparse data representation is preferred if there are implementations for both
    end
    
    properties (GetAccess = public)
        combinedPost
        noScore = -1;
        normalised = false;
        combinerInfo = '';
        nonDet = false;
    end
    
    properties (Abstract)
        label
        %shortLabel
    end

    methods
        function obj = AbstractCombiner(nAgents, K, targets, nClasses)
            
            if nargin > 3
                obj.nClasses = nClasses;
            end
            
            if exist('targets', 'var')            
                obj.setTargets(targets);
            end
            obj.nAgents = nAgents;
            obj.K = K;
        end
        
        function setTargets(obj, targets)
            obj.targets = targets-1;%transpose(targets);
        end
        
        function setProgressMonitor(obj, progressMonitor)
            obj.progressMonitor = progressMonitor;
        end
        
        function setId(obj, id)
            obj.id = id;
        end
        
        function id=getId(obj)
            id = obj.id;
        end
        %converts base outputs to the correct form or selects the correct
        %form if two are given. For the default version, any duplicate
        %scores are changed to the mean and the nAgents x nSamples matrix
        %is returned. 
        function baseData = correctBaseData(obj, baseOutputs)
            if iscell(baseOutputs) && (length(baseOutputs)==3 || length(baseOutputs)==4)

%                 [values, idxs] = unique([baseOutputs{1} baseOutputs{2}], 'rows', 'last');
%                 baseData = sparse(values(:,1), values(:,2), baseOutputs{3}(idxs));
                
                respCounts = sparse(baseOutputs{1}, baseOutputs{2}, 1);
                C3_avg = baseOutputs{3} ./ respCounts(sub2ind(size(respCounts),baseOutputs{1},baseOutputs{2}));
                C3_avg = round(C3_avg);
                obj.nAgents = max(baseOutputs{1});
                baseData = sparse(baseOutputs{1}, baseOutputs{2}, C3_avg);

                if length(baseOutputs)==4
                    obj.targets = baseOutputs{4};
                end
                
            else
                baseData = baseOutputs;                
            end
            
            obj.nAgents = size(baseData,1);
        end        
    end
    
    methods (Abstract)
        [combinedPost, agentRatings] = combineDecisions(obj, baseOutputs, clusters);
    end
end

