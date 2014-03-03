function [labelledTrainingData labelledTestData labels] = generateData(synthSet, expLabel, N, A)

synthSet = datageneration.expandSensorSettings(synthSet);

%GENERATE DATA
%create synthSet of training data
[D T] = datageneration.generateTwoClassData(synthSet);

labelledTrainingData = [D T];

dirName = datageneration.checkDataDir(synthSet);

dlmwrite(sprintf('%s%s_training_data.mat', dirName, expLabel), ...
    labelledTrainingData);

nInfSensors = synthSet.nInfSensors;
nNoninfSensors = synthSet.nNoninfSensors;

%expand deviations and means of informative sensors to all sensors
noninfDeviations = ones(2, 1) * synthSet.noninfDev;
deviations = [synthSet.deviations(:, 1:nInfSensors) ...
    10.*(ones(2,1)*rand(1,1))*ones(1,nNoninfSensors)];

noninfMeans = zeros(2, 1) + synthSet.noninfMean;
means = [synthSet.means(:, 1:nInfSensors) ...
    10.*(ones(2,1)*rand(1,1))*ones(1,nNoninfSensors)];

nSensors = synthSet.nSensors();

inf_sensors = zeros(1, nInfSensors);
noninf_sensors = zeros(1, nSensors - nInfSensors);
for i=1:nInfSensors
    inf_sensors(1, i) = i;
end
for i=nInfSensors+1:nSensors
    noninf_sensors(1, i-nInfSensors) = i;
end

%create a matrix for storing the test data and its labels and changepoint
%labels. Changepoints are represented by the last 2 columns; the first
%changepoint column is 0 or the number of a sensor that has been switched
%on; the second changepoint column is 0 or the number of the sensor that
%has just been switched off.

if ~exist('A', 'var')
    A = synthSet.nAgents;
end

labelledTestData = zeros(N, nSensors + 3);

%num corrupt labels
num_corrupt = round(rand(1, A) .* synthSet.pCorrupt_a) .* N;

%number of missing label intervals
%don't worry about overlaps for now.
num_missing = round(rand(1, A) .* synthSet.pMissing_a) .* N;

nSamples = N / synthSet.nDatasets;
if isempty(synthSet.defNoiseChanges)
    switchPoints = zeros(A, synthSet.switchPerAgent); %point at which an agent finishes switching between useless and informative
    switchEnds = zeros(A, synthSet.switchPerAgent); %point at which an agent stops changing
    period = floor(nSamples/(synthSet.switchPerAgent+1));
    for a=randperm(A)
        for s=1:synthSet.switchPerAgent
             %pick a data point to switch at. Exclude the first period of switch duration
            offset = 0;%round(a/A*period) - round(period/2);
            switchPoints(a,s) = round(s*period) - offset;
            switchEnds(a,s) = switchPoints(a,s) + synthSet.switchDuration;
        end
    end
    
    switchPoints = repmat(switchPoints, 1, synthSet.nDatasets);
    switchEnds = repmat(switchEnds , 1, synthSet.nDatasets);
    for d=1:synthSet.nDatasets
        switchPoints(:,(d-1)*synthSet.switchPerAgent+1:d*synthSet.switchPerAgent) = ...
            switchPoints(:,(d-1)*synthSet.switchPerAgent+1:d*synthSet.switchPerAgent) + (nSamples*(d-1));

        switchEnds(:,(d-1)*synthSet.switchPerAgent+1:d*synthSet.switchPerAgent) = ...
            switchEnds(:,(d-1)*synthSet.switchPerAgent+1:d*synthSet.switchPerAgent) + (nSamples*(d-1));   
    end
else
    switchPoints = synthSet.defNoiseChanges';
    switchEnds = synthSet.defNoiseChanges' + synthSet.switchDuration;
    
    if size(switchPoints,2)==1
        switchPoints = repmat(switchPoints, 1, synthSet.nDatasets);
        switchEnds = repmat(switchEnds , 1, synthSet.nDatasets);
    end
    for d=2:synthSet.nDatasets
        switchPoints(:,d:end) = switchPoints(:,d:end) + nSamples;
        switchEnds(:,d:end) = switchEnds(:,d:end) + nSamples;
    end    
    
    if size(switchPoints,1)<A
        switchPoints = [switchPoints; zeros(A-size(switchPoints,1),synthSet.nDatasets)];
    elseif size(switchPoints, 1)>A
        switchPoints = switchPoints(1:A, :);        
    end
    if size(switchEnds,1)<A
        switchEnds = [switchEnds; zeros(A-size(switchEnds,1),synthSet.nDatasets)];
    elseif size(switchEnds,1)>A
        switchEnds = switchEnds(1:A,:);
    end  

end



currentChanges = zeros(1, A);

pMode = synthSet.pMode;
pModeVariation = 0.5 - min(pMode);
minMode = min(pMode);
nFaultless = sum(pMode==minMode);

for n=1:N
    
    changepoint = [0, 0];
    
    %see if we should alter the informative sensor
    switch_inf = rand();
     %don't change near the start
    if (switch_inf < synthSet.pSwitch && n > 50) || sum(round(synthSet.definiteSwitch.*N)==n)>0
       display(sprintf('switching at %d', n));

       %randomly choose a sensor to switch off
       off_idx = ceil(rand() * nInfSensors);
       sensor_off = inf_sensors(1, off_idx);

       %choose one to switch on
       on_idx = ceil(rand() * (nSensors-nInfSensors) );
       sensor_on = noninf_sensors(1, on_idx);

       %switch new sensor to be the same as the old sensor
       means(:, sensor_on) = means(:, sensor_off);
       deviations(:, sensor_on) = deviations(:, sensor_off);   

       %switch the old one off
       means(:, sensor_off) = synthSet.noninfMean;
       deviations(:, sensor_off) = synthSet.noninfDev;

       %change the lists of sensors
       inf_sensors(1, off_idx) = sensor_on;
       noninf_sensors(1, on_idx) = sensor_off;
       
       changepoint = [sensor_on, sensor_off];
       
    else   
    %flip class labels? - don't do this if we have just switched sensors
        flip_inf_sensor = rand(1, nInfSensors);
        
        flip_idx = find(flip_inf_sensor < synthSet.pFlip_s);
        
        sensor_flip = inf_sensors(1, flip_idx);
            
        means([1 2], sensor_flip) = means([2 1], sensor_flip);
        deviations([1 2], sensor_flip) = means([2 1], sensor_flip);
        
        %do we need to record a changepoint?
    end
    
    if mod(n,nSamples)==0
        currentChanges(:) = 0;
        pMode = synthSet.pMode;
    end
        
   
    for a=1:nInfSensors
        changeNow = any(switchPoints(a,:)==n);
        stopChanging = any(switchEnds(a,:)==n);
        if changeNow
            if pMode(a) <= minMode %<= synthSet.pMode(a)
                currentChanges(a) = 1; %get stupider
                display([num2str(a) ' getting stupider']);
                
%                 if sum(pMode<=minMode)<nFaultless
%                     infa = find(pMode<=minMode);
%                     display([num2str(infa(1)) ' getting more informative*']);
%                     currentChanges(infa(1)) = -1;
%                 end
                
%                 switchPoints(a,:) = 0;
            elseif sum(pMode<=minMode)<nFaultless+1
                currentChanges(a) = -1; %get better
                display([num2str(a) ' getting more informative']);
                display(num2str(sum(pMode<=minMode)));
            else
                currentChanges(a) = -1; %get better
                display([num2str(a) ' getting more informative']);                
                %turn something else off
                infa = find(pMode<=minMode);
                currentChanges(infa(1)) = 1;
                display([num2str(infa(1)) ' getting stupider*']);
            end
        elseif stopChanging
            currentChanges(a) = 0;
        end
        
        pMode(a) = pMode(a) + (currentChanges(a)*pModeVariation) / synthSet.switchDuration;
        
        
        if pMode(a)>0.5
            pMode(a) = 0.5;
            currentChanges(a) = 0;
        elseif pMode(a) < minMode
            pMode(a) = minMode;
            currentChanges(a) = 0;
        end
        
    end
            
    [point classification] = datageneration.generateDataPoint(nSensors, ...
        nInfSensors, means, deviations, synthSet.p_c1, pMode, synthSet.class1ErrorOnly);
                
    labelledTestData(n, :) = [point classification changepoint];
end

labels_for_agents = labelledTestData(:, nSensors+1) * ones(1, A);
%corrupt some labels

%randomly choose the agents to corrupt
%A = ceil(rand(num_corrupt ,1)*A);

%corruptLabels has two columns. First column records the sample number
%where a string of corrupted labels starts. Second column contains a matrix
%with the agent number that was corrupted, with length corresponding to
%number of labels corrupted.
corruptLabels = zeros(sum(num_corrupt), 2);
for a=1:A
    
    if num_corrupt(a) == 0
        continue
    end
    
    corrupt_x = ceil( rand(num_corrupt(a), 1) * N );
    
    %is this correct? The labelledTestData has no range in the rows
    labels_for_agents(corrupt_x, a) = ones(num_corrupt,1) - ...
        labelledTestData(corrupt_x, nSensors+1);
    
    start_idx = sum(num_corrupt(1:a-1)) + 1;
    corruptLabels(start_idx:num_corrupt(a), 1) = corrupt_x;
    corruptLabels(start_idx:num_corrupt(a), 2) = ones(1,num_corrupt(a))*a;
end

synthSet.corruptLabels = corruptLabels;

expanded_pFlip_a = ones(N,1)*synthSet.pFlip_a;
flippedLabels = rand(N, A);

%flipped labels contains 0s where label is not flipped
flippedLabels( flippedLabels > 1-expanded_pFlip_a ) = 1;
flippedLabels( flippedLabels <= 1-expanded_pFlip_a) = 0;
for a=1:A
    flips = find(flippedLabels(:,a));
    for f=1:length(flips)
        labels_for_agents(flips(f):N, a) = ...
            1 - labels_for_agents(flips(f):N, a); 
    end
end

synthSet.flippedLabels = flippedLabels;

missingLabels = zeros(sum(num_missing)*2, 3);
for a=1:A
    
    if num_missing(a) == 0
        continue
    end
    
    missing_x = ceil( rand(num_missing(a), 1) .* N );
    missing_lengths = ceil( randn(num_missing(a), 1) * dev_ml_a(a) + mean_ml_a(a));
    missing_xend = missing_x + missing_lengths;
    %correct any end values that have gone over the end of test data
    missing_xend(find(missing_xend>N), 1) = N;
    
    start_idx = sum(num_missing(1:a-1))*2 + 1;
    missingLabels(start_idx:num_missing(a), 1) = missing_x;
    missingLabels(start_idx:num_missing(a), 2) = ones(1,num_missing(a));
    missingLabels(start_idx:num_missing(a), 3) = ones(1,num_missing(a)) * a;
    
    start_idx = sum(num_missing(1:a-1))*2 + num_missing(a) + 1;
    missingLabels(start_idx:num_missing(a), 1) = missing_xend;
    missingLabels(start_idx:num_missing(a), 2) = zeros(1,num_missing(a));
    missingLabels(start_idx:num_missing(a), 3) = ones(1,num_missing(a)) * a;

    for x=1:length(missing_x)
        labels_for_agents(missing_x(x, 1):missing_xend(x, 1), a) = -1;
    end    
end

synthSet.missingLabels = missingLabels;

dlmwrite(sprintf('%s%s_test_data.mat', synthSet.getDataDir(), expLabel), ...
    labelledTestData);

fprintf('Generated %d training data points and %d test data points\n', ...
    synthSet.nTrainingSamples, N);

labels = labelledTestData(1:N, nSensors+1);

end