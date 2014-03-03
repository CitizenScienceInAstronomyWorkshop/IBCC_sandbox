function synthSet = expandSensorSettings( synthSet )

%expandSensorSettings Fill in defaults where there are missing values

    datageneration.defaultSensorSettings;

    specifiedDevs = size(synthSet.deviations, 2);
    if specifiedDevs < synthSet.nInfSensors
        synthSet.deviations = [synthSet.deviations ...
            defDev * ones(1, synthSet.nInfSensors - specifiedDevs) ...
            ];
        
    elseif specifiedDevs > synthSet.nInfSensors
        synthSet.deviations = synthSet.deviations(:, 1:synthSet.nInfSensors);
    end

    specifiedsynthSet.means = size(synthSet.means, 2);
    if specifiedsynthSet.means < synthSet.nInfSensors
        
        synthSet.means = [synthSet.means ...
            defMean*ones(1, synthSet.nInfSensors - specifiedsynthSet.means)...
            ];
        
    elseif specifiedsynthSet.means > synthSet.nInfSensors
        synthSet.means = synthSet.means(:, 1:synthSet.nInfSensors);
    end

    specifiedMissing = size(synthSet.pMissing_a, 2);
    if specifiedMissing < synthSet.nAgents
        
        synthSet.pMissing_a = [synthSet.pMissing_a ...
            ones(1, synthSet.nAgents-specifiedMissing).*synthSet.pMissing];
        
    elseif specifiedMissing > synthSet.nAgents
        synthSet.pMissing_a = synthSet.pMissing_a(:, 1:synthSet.nAgents);
    end

    specifiedMissing = size(synthSet.mean_ml_a, 2);
    if specifiedMissing < synthSet.nAgents;
        
        synthSet.mean_ml_a = [synthSet.mean_ml_a ...
            ones(1, synthSet.nAgents-specifiedMissing).*synthSet.meanMissingLength];
        
    elseif specifiedMissing > synthSet.nAgents
        synthSet.mean_ml_a = synthSet.mean_ml_a(:, 1:synthSet.nAgents);
    end

    specifiedMissing = size(synthSet.dev_ml_a, 2);
    if specifiedMissing < synthSet.nAgents;
        
        synthSet.dev_ml_a = [synthSet.dev_ml_a ...
            ones(1, synthSet.nAgents-specifiedMissing).*synthSet.devMissingLength];
        
    elseif specifiedMissing > synthSet.nAgents
        synthSet.dev_ml_a = synthSet.dev_ml_a(:, 1:synthSet.nAgents);
    end

    specified_corrupt = size(synthSet.pCorrupt_a, 2);
    if specified_corrupt < synthSet.nAgents
        
        synthSet.pCorrupt_a = [synthSet.pCorrupt_a ...
            ones(1, synthSet.nAgents-specified_corrupt).*synthSet.pCorrupt];
        
    elseif specified_corrupt > synthSet.nAgents
        synthSet.pCorrupt_a = synthSet.pCorrupt_a(:, 1:synthSet.nAgents);
    end

    specified_flip = size(synthSet.pFlip_a, 2);
    if specified_flip < synthSet.nAgents
        
        synthSet.pFlip_a = [synthSet.pFlip_a ...
            ones(1, synthSet.nAgents-specified_flip).*synthSet.pFlipLabel];
        
    elseif specified_flip > synthSet.nAgents
        synthSet.pFlip_a = synthSet.pFlip_a(:, 1:synthSet.nAgents);
    end
    
    specified_flip = size(synthSet.pFlip_s, 2);
    if specified_flip < synthSet.nInfSensors
        
        synthSet.pFlip_s = [synthSet.pFlip_s ...
            ones(1, synthSet.nInfSensors-specified_flip).*synthSet.pFlipSensor];
        
    elseif specified_flip > synthSet.nInfSensors
        synthSet.pFlip_s = synthSet.pFlip_s(:, 1:synthSet.nInfSensors);
    end
    
    

end

