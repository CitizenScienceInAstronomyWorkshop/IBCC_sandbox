if ~exist('rootDir','var')
    rootDir = [];
end
expSet = settings.ExpSettings(rootDir, 'gzsnData.csv');

expSet.nRepeats = 1;
expSet.topSaveDir = '/homes/49/edwin/matlab/data/cyber/combiner_output';%exp2its';%
expSet.expLabel = 'bcc';
expSet.outputDir = '/homes/49/edwin/matlab/data/cyber/combiner_output';
expSet.saveAs = '';

expSet.voteThreshold = 0.5;
expSet.nScores = 2;
expSet.nClasses = 8;

bccSet = settings.BccSettings();
bccSet.Alpha = [0.5 0.5; 0.5 0.5; 0.5 0.5; 0.5 0.5; 0.5 0.5; 0.5 0.5; 0.5 0.5; 0.5 0.5];
bccSet.nu = {[1000 1000 1000 1000 1000 1000 1000 1000]};

expSet.noScore = -1;
bccSet.maxScore = 1;
bccSet.minScore = 0;

bccSet.debug = true;

bccSet.scoreMap = [];
