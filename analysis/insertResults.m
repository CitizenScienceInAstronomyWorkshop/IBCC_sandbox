clear
import settings.*;

%experiment 1
bccSettings;
load(sprintf('%s%s', '/homes/49/edwin/matlab/combination/data/bcc3/comparison/', expSettings.multCombTestFile));

%experiment 2
%bccSettingsHard;
%load(sprintf('%s%s', '/homes/49/edwin/matlab/combination/data/bcc4HardWM/comparisonHard/', expSettings.multCombTestFile));

%experiment 3
%bccSettingsExp3;
%load(sprintf('%s%s', '/homes/49/edwin/matlab/combination/data/bccExp3WM/comparisonMinorityInformed/', expSettings.multCombTestFile));

%exp 4
%bccSettingsInversions;
%load(sprintf('%s%s', '/homes/49/edwin/matlab/combination/data/bccInvWM/comparison/', expSettings.multCombTestFile));

newResults = combinedPostKnowns;

load(sprintf('%s%s', expSettings.getCombinerDir(), expSettings.multCombTestFile));

%insert into first setting only

for d=1:10    
    for i=1:length(expSettings.propKnown)
        oldResults = combinedPostKnowns{d, i};
        resToInsert = newResults{d, i};
        
        resToInsert = round(resToInsert);
        
        results = zeros(size(oldResults, 1), size(oldResults, 2), 5);
        results(:, :, 1:4) = oldResults;
        results(:, :, 5) = resToInsert(:, :, 4);
    
        combinedPostKnowns{d, i} = results;
    end
end

save(sprintf('%supdatedNB_%s', expSettings.getCombinerDir(), expSettings.multCombTestFile), 'combinedPostKnowns');