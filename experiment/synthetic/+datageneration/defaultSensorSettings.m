
defMean = [1; -1]; %[-3 3];
defDev = [2; -2]; %[1 1];

pMissing = 0;
meanMissingLength = 10; %period of "I don't know" data points
devMissingLength = 3;

%switch the label that a particular agent sees for a single data point
%how useful is this?
pCorrupt = 0;

%switch all the labels beyond this point for a particular agent
pFlipLabel = 0; %0.005;

%flip the sensor permanently
pFlipSensor = 0;