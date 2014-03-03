classdef GzRunner < ExpRunner
    %GZRUNNER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = GzRunner(expSettings, baseScores, trainingLabels, allLabels)
            if isempty(allLabels)
                expSettings.nSamples = max(baseScores{2});
            else
                expSettings.nSamples = length(allLabels);
            end
            
            obj@ExpRunner(expSettings);
            
            obj.loadData = false;
            if ~isempty(allLabels)
                if size(allLabels, 1) > 1 && size(allLabels, 2)==1
                    allLabels = allLabels';
                end
                obj.labels = allLabels;
                expSettings.nSamples = length(obj.labels);
            end
            
            if ~isempty(trainingLabels)
                obj.trainingLabels = trainingLabels;
            end
            
            if ~isempty(baseScores)
                obj.basePost = [];
                obj.baseScoreSet = baseScores;
            end
            obj.postScores = false;
        end        
    end
    
end

