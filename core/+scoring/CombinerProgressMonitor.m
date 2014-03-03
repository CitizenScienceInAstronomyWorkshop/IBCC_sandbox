classdef CombinerProgressMonitor < handle
    %COMBINERPROGRESSMONITOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        labels
        testLabels
        testIdxs
        error
        records
        nRepeats = 0;
        nKnown
        
        method = 'AUC';%'absolute error';
        
        prevMax = inf;
    end
    
    methods
        function obj = CombinerProgressMonitor(labels, nCombiners)
            obj.labels = labels;
            obj.error = cell(nCombiners, 1);
            obj.records = cell(nCombiners, 1);
            obj.nRepeats = zeros(1,nCombiners);
        end
        
        function updateLabels(obj, labels, testIdxs, nKnown)
            obj.labels = labels;
            if nargin > 2
                obj.testLabels = zeros(1, length(obj.labels));
                obj.testLabels(testIdxs) = obj.labels(testIdxs);
                obj.testIdxs = testIdxs;
                obj.nKnown = nKnown;
            end
        end
        
        function combinerIteration(obj, currentResults, combinerId, nIt, nKnown)
            %nKnown assumes the nKnown targets are at the start.
            if nargin < 5
                nKnown = 0;
            end
            
            if size(currentResults,1)>1
                currentResults = currentResults(2:end,:);
            end
            
            if strcmp(obj.method, 'pass value in')
                %do nothing
                testLabels = [];
            elseif ~isempty(obj.nKnown) && ~isempty(obj.testLabels)
                nKnown = obj.nKnown;
                testResults = currentResults(:,obj.testIdxs);
                testLabels = obj.testLabels(obj.testIdxs);
                nUnknown = sum(obj.testLabels==0);
            else
                nUnknown = length(currentResults) - nKnown;            
                testResults = currentResults(nKnown+1:end);
                testLabels = obj.labels(nKnown+1:length(currentResults));
            end
                  
            if size(testLabels,1) > size(testLabels, 2)
                testLabels = testLabels';
            end
            
            if strcmp(obj.method, 'absolute error')           
                currentError = sum(abs(testResults-testLabels));
                        
                if  ~isempty(obj.nKnown) && ~isempty(obj.testLabels)
                    if currentError > 0.5 * (nKnown+nUnknown)
                        currentError = (nKnown+nUnknown)-currentError;
                    end                 
                else
                    currentError = currentError * length(obj.labels) / nUnknown;
                    if currentError > 0.5 * (nKnown + nUnknown)
                        currentError = length(obj.labels)-currentError;
                    end               
                end   
            elseif strcmp(obj.method, 'mean squared error')           
                currentError = sum((testResults-testLabels).^2) ./ length(testLabels);
            elseif strcmp(obj.method, 'cross entropy')
                testLabels = testLabels - 1;
                testResults = testResults -1;
                
                
                currentError = -testLabels.*log(testResults) -(1-testLabels).*log(1-testResults);
                
                
                
%                 currentError(testResults==0&testLabels==0) = 0;
%                 currentError(testResults==1&testLabels==1) = 0;
%                 currentError(testResults==0&testLabels==1) = Inf;
%                 currentError(testResults==1&testLabels==0) = Inf;
                currentError = sum(currentError);
            elseif strcmp(obj.method, 'AUC')
                nClasses = size(testResults,1);
                currentError = zeros(nClasses+1, 1);
                noGraphs  = true;
                for j=1:nClasses
                    currentError(j) = graphs.ClassifierPerformanceGraph.drawRoc(testResults(j,:), ...
                        testLabels-1, {'topic'}, false, noGraphs, false);
                end
                %also include the entropy
                
                logTestResults = log(testResults);
                logTestResults(isnan(logTestResults) | isinf(logTestResults)) = 0;
                
                currentError(nClasses+1) = -sum(sum(testResults .* logTestResults));
                
            elseif strcmp(obj.method, 'pass value in')
                currentError = currentResults;
            end
            
            if nIt<=obj.prevMax
                newRep = true;
            else
                newRep = false;
            end     
            obj.prevMax = nIt;
            
            if nargin < 3 || nIt > size(obj.error{combinerId},2)
                if newRep
                    obj.nRepeats(combinerId) = obj.nRepeats(combinerId) + 1;
                end
                if newRep
                    nRepIts = nIt;
                else
                    nRepIts = nIt - size(obj.error{combinerId},2);
                end
                currentError = repmat(currentError, 1, nRepIts);
                
                obj.error{combinerId} = [obj.error{combinerId} currentError];
                if obj.nRepeats(combinerId)==0
                    obj.nRepeats(combinerId) = 1;
                end
            else
                if newRep
                    obj.nRepeats(combinerId) = obj.nRepeats(combinerId) + 1;
                end    
                
                if newRep
                    nRepIts = nIt;
                else
                    nRepIts = 1;
                end
                currentError = repmat(currentError, 1, nRepIts);
                
                obj.error{combinerId}(:,nIt-nRepIts+1:nIt) = obj.error{combinerId}(:,nIt-nRepIts+1:nIt) .* ...
                    (obj.nRepeats(combinerId)-1) + currentError;
                obj.error{combinerId}(:,nIt-nRepIts+1:nIt) = obj.error{combinerId}(:,nIt-nRepIts+1:nIt) ./ ...
                    obj.nRepeats(combinerId);
            end
            if size(obj.records,2)<obj.nRepeats(combinerId)
                obj.records{combinerId,obj.nRepeats(combinerId)} = currentError;
            else
                obj.records{combinerId,obj.nRepeats(combinerId)} = [obj.records{combinerId,obj.nRepeats(combinerId)} currentError];
            end
        end
        
        function SD = calculateDeviation(obj)
            %calculates the standard deviation at each iteration
            
            SD = cell(length(obj.error),1);
            
            for c=1:size(obj.records,1)
                deviations = zeros(size(obj.error{c}));
                for rep=1:obj.nRepeats(c)
                    deviations = deviations + (obj.records{c,rep} - obj.error{c}).^2;
                end
                SD{c} = (deviations ./ obj.nRepeats(c)).^0.5;
            end
        end
    end
    
end

