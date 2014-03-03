function [snBaseOutputs, snRawData, labels, typeLabels, typeAssets] = ...
    reloadGZSNData(loadDataFromFile, dataIsAnnotations, filename, snRawData)

if loadDataFromFile && isempty(snRawData)
    if strcmp(filename(end-3:end), '.mat')
        load(filename);
        snRawData = snData;
    else %assume .csv format
        fid = fopen(filename);
        snRawData = textscan(fid, '%d %d %d %d %d', 'delimiter', ',');
        fclose(fid);
        snData = snRawData;
    end
else
    snData = snRawData;
end

%see extract_useful_data.sql
% cells of snData
%1.  Task ID (MWO)
%2. Agent ID (Site)
%3. Score (Choice)
%4. Label (Bench Result)
%5. Symptom

%use class labels as our label set. Ignore uncertain labels.
correctScoreIdxs = find(snData{4}==1);
correctScores = snData{3}(correctScoreIdxs);
[labelIds, labelIdxs] = unique(snData{1}(correctScoreIdxs));
labels = correctScores(labelIdxs) + 1;

snBaseOutputs = cell(3,1);
snBaseOutputs{1} = double(snData{2} + 1); %agents
snBaseOutputs{2} = double(snData{1} + 1); %assets
snBaseOutputs{3} = double(snData{3} + 1); %scores

typeAssets = [];
typeLabels = [];
end