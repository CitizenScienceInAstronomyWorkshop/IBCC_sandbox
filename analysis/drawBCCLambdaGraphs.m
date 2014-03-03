clear

import settings.*;

bccLambdaSettings;

load(sprintf('%sibccGibbsTests-propKnown-nSamples100.mat', expSettings.getCombinerDir()));
labelledTestData = dlmread(sprintf('%sdlrCombiners_test_data.mat', expSettings.getDataDir()));
labels = labelledTestData(:, expSettings.nSensors()+1);
labels = reshape(labels, expSettings.nSamples, expSettings.nDatasets)';

varSettings = expSettings.propKnown;

graphMaker = graphs.VariableSettingGraph('Proportion of known targets', expSettings.nSamples);

[resKnownsMeans, resKnownsDev, resKnownsMin, resKnownsMax] = ...
    graphMaker.error2dRuns(combinedPostKnowns, labels, true);
graphMaker.drawErrorGraphDev(resKnownsMeans, resKnownsDev, varSettings, [], false);
graphMaker.drawErrorGraphMinMax(resKnownsMeans, resKnownsMin, resKnownsMax, ...
    varSettings, [], false);

load(sprintf('%sibccGibbsTests2-lambda-nSamples100.mat', expSettings.getCombinerDir()));

[resLambdaSymMeans, resLambdaSymDev, resLambdaSymMin, resLambdaSymMax, ...
    resLambdaMagMeans, resLambdaMagDev, resLambdaMagMin, resLambdaMagMax] = ...
    graphMaker.error3dRuns(combinedPostLambda, labels, true);

%draw graphs to show affect of changing lambda symmetry. Plot each
%magnitude value of lambda separately.
resLambdaSymMeans = sum(resLambdaSymMeans, 1) ./ size(resLambdaSymMeans, 1);
resLambdaSymDev = sum(resLambdaSymDev, 1) ./ size(resLambdaSymDev, 1);
resLambdaSymMin = min(resLambdaSymMin, [], 1);
resLambdaSymMax = max(resLambdaSymMax, [], 1);

graphMaker.label = 'totals for symmetry between diags and non-diags in lambda';
graphMaker.drawErrorGraphDev(resLambdaSymMeans, resLambdaSymDev, ...
    expSettings.lambdaSym, [], false);
graphMaker.drawErrorGraphMinMax(resLambdaSymMeans, resLambdaSymMin, resLambdaSymMax, ...
    expSettings.lambdaSym, [], false);

%changing lambda magnitude - plot different symmetry levels separately
resLambdaMagMeans = sum(resLambdaMagMeans, 1) ./ size(resLambdaMagMeans, 1);
resLambdaMagDev = sum(resLambdaMagDev, 1) ./ size(resLambdaMagDev, 1);
resLambdaMagMin = min(resLambdaMagMin, [], 1);
resLambdaMagMax = max(resLambdaMagMax, [], 1);

graphMaker.label = 'magnitude of lambda - each plot corresponds to different symmetry level';
graphMaker.drawErrorGraphDev(resLambdaMagMeans, resLambdaMagDev, ...
    expSettings.lambdaMag, expSettings.lambdaSym, false);
graphMaker.drawErrorGraphMinMax(resLambdaMagMeans, resLambdaMagMin, resLambdaMagMax,...
    expSettings.lambdaMag, expSettings.lambdaSym, false);