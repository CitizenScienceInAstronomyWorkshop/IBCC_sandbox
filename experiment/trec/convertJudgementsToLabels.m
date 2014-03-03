fileMapFile = '/homes/49/edwin/data/trec/finalOutput/MatrixTopic2000/Topic-Matrix/fileMap.txt';
judgeFile = '/homes/47/parg/nsdata/Data/classifier_combination/trec/trec-2012-trat-adjudicated-judgments-Oct-11-2012';

fid = fopen(fileMapFile);
fileMapCells = textscan(fid, '%d %s %s %s', 'Delimiter', '/_., ');

fid = fopen(judgeFile);
judgeCells = textscan(fid, '%d %d %s %d', 'Delimiter', ' ');


fileIdxs = fileMapCells{1} + 1; %starts at zero
fileMapNames = fileMapCells{3};

posIdxs = judgeCells{4}~=0;
judgeLabels = zeros(length(judgeCells{1}), 1);
judgeLabels(posIdxs) = judgeCells{1}(posIdxs);
judgeNames = judgeCells{3};

labels = zeros(length(fileIdxs),1);

clashIdxs = [];

for i=1:length(judgeNames)
    
    if judgeLabels(i)==0
        continue;
    end
    
    name = judgeNames{i};
    fileIdx = fileIdxs(strcmp(name, fileMapNames));
    
    if labels(fileIdx)==0
        labels(fileIdx) = judgeLabels(i);
    else
        display(['clash: ' num2str(labels(fileIdx)) ', ' num2str(judgeLabels(i))]);
        clashIdxs = [clashIdxs fileIdx];
    end
end

save('/homes/49/edwin/data/trec/judgementsAsArray.mat', 'labels');

[uniqueLabels, ~, mappedLabels] = unique(labels);

save('/homes/49/edwin/data/trec/judgementsAsMappedArray.mat', 'mappedLabels');

mappedLabelsSkipClashes = mappedLabels;
mappedLabelsSkipClashes(clashIdxs) = 0;

save('/homes/49/edwin/data/trec/judgementsAsMappedArray_skipClashes.mat', 'mappedLabelsSkipClashes');
