function mergeResults( expSettings, secondDataDir, mergedDataDir, exMethodIdx1, exMethodIdx2 )
%exMethodIdx - indexes of combination method to exclude.

[labels1, results1, combMethods1] = loadMultTestResults(expSettings);
expSettings.topSaveDir = secondDataDir;%'/homes/49/edwin/matlab/combination/data/bcc3prs';
[labels2, results2, combMethods2] = loadMultTestResults(expSettings);

combinedPostKnowns = cell(expSettings.nDatasets, expSettings.nRepeats);

nCombMethods1 = length(combMethods1);
nCombMethods2 = length(combMethods2);

methodIdx1 = 1:nCombMethods1;
methodIdx2 = 1:nCombMethods2;

for d=1:expSettings.nDatasets
    for i=1:length(expSettings.propKnown)
        
        cleanResults1 = results1{d,i};
        cleanResults1 = cleanResults1(:, :, methodIdx1( ~ismember(methodIdx1, exMethodIdx1)) );
        if size(cleanResults1, 1) == 1 && expSettings.nRepeats > 1
            display(sprintf('%s: Only found one repetition - assuming these methods are deterministic so each results set can be duplicated as they are the same for every repeat', expSettings.getDataDir()));
            expResults1 = zeros(expSettings.nRepeats, expSettings.nSamples, size(cleanResults1, 3));
            for r=1:expSettings.nRepeats
                expResults1(r, :, :) = cleanResults1;
            end
            cleanResults1 = expResults1;
        end
        
        cleanResults2 = results2{d,i};
        cleanResults2 = cleanResults2(:, :, methodIdx2( ~ismember(methodIdx2, exMethodIdx2)) );
        if size(cleanResults2, 1) == 1 && expSettings.nRepeats > 1
            display(sprintf('%s: Only found one repetition - assuming these methods are deterministic so each results set can be duplicated as they are the same for every repeat', secondDataDir));
            expResults2 = zeros(expSettings.nRepeats, expSettings.nSamples, size(cleanResults2, 3));
            for r=1:expSettings.nRepeats
                expResults2(r, :, :) = cleanResults2;
            end
            cleanResults2 = expResults2;
        end
        
        
        combinedPostKnowns{d, i} = zeros(expSettings.nRepeats, ...
            size(cleanResults1, 2), size(cleanResults1, 3)+size(cleanResults2, 3));
        
        combinedPostKnowns{d, i}(:, :, 1:size(cleanResults1, 3)) = cleanResults1; 
        combinedPostKnowns{d, i}(:, :, size(cleanResults1, 3)+1:end) = cleanResults2;
    end
end

combMethods1 = combMethods1( ~ismember(methodIdx1, exMethodIdx1) );
combMethods2 = combMethods2( ~ismember(methodIdx2, exMethodIdx2) );

combMethods = cell(1, length(combMethods1)+length(combMethods2));
combMethods(1:length(combMethods1)) = combMethods1;
combMethods(length(combMethods1)+1:length(combMethods1)+length(combMethods2)) = combMethods2;

if ~exist(mergedDataDir, 'file')
    mkdir(mergedDataDir);
end
if ~exist(sprintf('%s/%s', mergedDataDir, expSettings.expLabel), 'file')
    mkdir(sprintf('%s/%s', mergedDataDir, expSettings.expLabel));
end
if ~exist(sprintf('%s/baseData', mergedDataDir), 'file')
    mkdir(sprintf('%s/baseData', mergedDataDir));
end

save(sprintf('%s/%s/combMethods-%s', mergedDataDir, expSettings.expLabel, expSettings.multCombTestFile), ...
    'combMethods');
save(sprintf('%s/%s/%s', mergedDataDir, expSettings.expLabel, expSettings.multCombTestFile), ...
    'combinedPostKnowns');
labelsFile = sprintf('%s/baseData/%s_labels.mat', mergedDataDir, expSettings.expLabel);
dlmwrite(labelsFile, labels1');
end

