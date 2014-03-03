function [ newLabels, responders ] = getSimResponses(agents, workersToRequest, ...
    samplesToRequest, excludedWorkers, qRels)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    %Get the requested labels from the prepared set
    newLabels = [];
    responders = [];
    
    nWorkers = length(agents);

    for t=1:length(samplesToRequest)
        w = workersToRequest(t);
        i = samplesToRequest(t);
        if w==0
            while w==0 || ismember(w, excludedWorkers)
                w = randi(nWorkers,1,1);
            end
            resp = agents{w}.getResponse(qRels,i);
            if resp~=-1
                newLabels = [newLabels; w i resp 0 0];
%                     workerQueue = [workerQueue newWorker];
            end
        else
            resp = agents{w}.getResponse(qRels,i);
            if resp~=-1
                newLabels = [newLabels; w i resp 0 0];
%                     workerQueue = [workerQueue w];
            end
        end
        responders = [responders; w];
    end
end

