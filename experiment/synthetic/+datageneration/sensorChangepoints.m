function changes = sensorChangepoints(synthSet, nSamples, labelledTestData)
    changepoints_plus = [];
    changepoints_minus = [];
    changepoints_x = [];

    %record initial states of sensors
    num_changepoints = 0;
    for i=1:synthSet.nInfSensors
        num_changepoints = num_changepoints + 1;
        changepoints_x(1, num_changepoints) = 1;
        changepoints_plus(1, num_changepoints) = i;
    end
    
    nSensors = synthSet.nSensors();

    %record subsequent changes to sensors
    changes_start = num_changepoints+1;
    for i=1:nSamples

        if labelledTestData(i, nSensors+2) > 0
            num_changepoints = num_changepoints + 1;

            changepoints_plus(1, num_changepoints) ...
                = labelledTestData(i, nSensors+2);

            changepoints_minus(1, num_changepoints) ...
                = labelledTestData(i, nSensors+3);

            changepoints_x(1, num_changepoints) = i;
        end
        
    end
    
    changes = {changes_start changepoints_plus changepoints_minus changepoints_x num_changepoints};
    
    %NEW VERSION
    %don't need to return num changepoints or changes start
    
    sparseChangesPlus = transpose(labelledTestData(:, nSensors+2));
    sparseChangesMinus = transpose(labelledTestData(:, nSensors+3));
    
    changesPlusN = sparseChangesPlus>0;
    changesMinusN = sparseChangesMinus>0;
    
    %append the initial states as if they were sensors becoming active
    changesPlusN = [ones(1, synthSet.nInfSensors) find(changesPlusN)];
    changesMinusN = find(changesMinusN);
      
    
    changepointsPlus = sparseChangesPlus(changesPlusN);
    %add in the sensor IDs for the initial informative sensors
    changepointsPlus(1:synthSet.nInfSensors) = 1:synthSet.nInfSensors;
    
    changepointsMinus = sparseChangesMinus(changesMinusN);
    
    changes = {changesPlusN changepointsPlus changesMinusN changepointsMinus};
end