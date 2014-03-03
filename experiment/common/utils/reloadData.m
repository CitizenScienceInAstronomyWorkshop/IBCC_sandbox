function [snBaseOutputs, snRawData, labels, assetIds] = ...
    reloadData(loadDataFromFile, dataIsAnnotations, filename, snRawData)

if loadDataFromFile && isempty(snRawData)
    display(['Reading input data from file ' filename]);

    fid = fopen(filename);
    snRawData = textscan(fid, '%d %d %d %d', 'delimiter', ',');
    fclose(fid);
    snData = snRawData;
else
    snData = snRawData;
end

%see extract_useful_data.sql
% cells of snData:
%1. user ID
%2. asset ID
%3. class
%4. score

%If we are using a full list of annotations use this version
if nargin > 1 && dataIsAnnotations
    %classifications only, not the individual annotations
    classData = snData;

    [classData{2}, uniqueIdxs] = unique(snData{2});
    for i=1:length(classData)
        if i==2
            continue %done it already
        end
        if ~isempty(snData{i})
            classData{i} = snData{i}(uniqueIdxs);
        end
    end

    scores = classData{4};
    assets = classData{2};
    agents = classData{1};
    classes = classData{3};
end

%use class labels as our label set. Ignore uncertain labels.
[assetIds uniqueAssetIdxs assetIdxs] = unique(assets);
[vals uniqueAgentIdxs agentIdxs] = unique(agents); %get the unique asset IDs - previously they were duplicated across each zoo classification
labels = classes(uniqueAssetIdxs);

%This is probably incorrect if the type labels are no wholly reliable -
%include these as reliable agents only
% %Currently, treat those with no from type labelling as definite no only
% negTypes = assetTypeLabels==-1 & assetClassLabels==0;
% labels(negTypes) = 1; 

snBaseOutputs = cell(3,1);
snBaseOutputs{1} = double(agentIdxs); %agents
snBaseOutputs{2} = double(assetIdxs); %assets
snBaseOutputs{3} = double(scores); %scores
end