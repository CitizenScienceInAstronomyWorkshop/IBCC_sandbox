classdef BaeSimRunner < ExpRunner
    %GZRUNNER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = BaeSimRunner(expSettings, baseScores, trainingLabels, allLabels, Alpha)
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
            obj.noScore = 0;
            obj.maxScore = 10;
            obj.minScore = 1;
            obj.voteThreshold = 0;
            obj.nScores = 10; 
            obj.nClasses = 10;
            
            obj.targetsAsAgents = [];%[1 2];
            if nargin < 5 || isempty(Alpha)

              obj.Alpha0 = ones(obj.nClasses, obj.nScores) .* 0.7 ./ obj.nScores;
              obj.Alpha0(1:obj.nClasses, 1:obj.nClasses) = 0.3;
            else
                obj.Alpha0 = Alpha;
            end
            obj.IbccMapFunction = @dummyMap;

        end        
    end
    
end

