function [snBaseOutputs, snRawData, labels, typeLabels, typeAssets, assetIds] = ...
    reloadGZSNData(loadDataFromFile, dataIsAnnotations, filename, snRawData, interpretNoTypes)

if nargin < 5
    interpretNoTypes = true;
end

if loadDataFromFile && isempty(snRawData)
    if strcmp(filename(end-3:end), '.mat')
        load(filename);
        snRawData = snData;
    else %assume .csv format
        
        display(['Reading input data from file ' filename]);
        
        fid = fopen(filename);
        snRawData = textscan(fid, '%d %d %d %s %s %d', 'delimiter', ',');
        fclose(fid);
        snData = snRawData;
    end
else
    snData = snRawData;
end

%see extract_useful_data.sql
% cells of snData:
%1. annotation ID
%2. classification ID
%3. user ID
%4. asset ID
%5. asset name
%6. PTF ID
%7. PTF type value 
%8. PTF class value
%9. total score
%10. no. classifications for this asset
%11. answer ID
%12. task ID
%13. answer score
%14. updated annotation at...
if length(snData)<14 %fill in some empty slots, assume 1, 5 and 6 have been missed
    snData = [{[]} snData(1:3) {[] []} snData(4:end)];
end    

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

    scores = classData{9};
    assets = classData{4};
    agents = classData{3};
    ptfType = classData{7};
    ptfClass = classData{8};
else
    scores = snData{7};
    assets = snData{4};
    agents = snData{2};
    ptfType = snData{6};
    ptfClass = snData{7};
end

%remove crappy data where score is impossible
validIdx = find(scores==-1 | scores==1 | scores==3);
scores = scores(validIdx);
assets = assets(validIdx);
agents = agents(validIdx);
ptfType = ptfType(validIdx);
ptfClass = ptfClass(validIdx);

%use class labels as our label set. Ignore uncertain labels.
[assetIds uniqueAssetIdxs assetIdxs] = unique(assets);
[vals uniqueAgentIdxs agentIdxs] = unique(agents); %get the unique asset IDs - previously they were duplicated across each zoo classification
assetTypeLabels = ptfType(uniqueAssetIdxs);
assetClassLabels = ptfClass(uniqueAssetIdxs);

%Translate the PTF labels into scores - only -1, 1 and 3 are allowed.
classCounts = sparse(length(assetIds), 3); %record the number of PTF entries supporting each classification
typeCounts = sparse(length(assetIds), 2);
numTypeData = sparse(length(assetIds), 1); %count the number of PTF entries saying an asset is a transient
numClassData = sparse(length(assetIds), 1); %count the number of PTF entries saying an asset is a SN
nClassOnly = 0; %number for which we have no PTF type 

for i=1:length(assetIds)
    type = assetTypeLabels(i);
    type = type{1};
    tokens = textscan(type, '%s');
    if strcmp(type, '""')
        tokens = {[]};
    end
    for t=1:length(tokens{1})
        token = tokens{1}{t};
        
        if strcmp(token, '"')
            continue
        end
        
        if ~isempty(strfind(token, 'Transient'))
            typeCounts(i, 2) = typeCounts(i, 2) + 1;
        else
%             display(['token: ' token])
            typeCounts(i, 1) = typeCounts(i, 1) + 1;
        end
    end
    [maxScore, maxIdx] = max(typeCounts(i,:));
    maxFlags = typeCounts(i, :) == maxScore;
    % Non-transient -> -1, transient -> 1, >1 transient labellings -> 3
    if sum(maxFlags) == 1
        if maxIdx==1
            numTypeData(i) = -1;
        elseif maxScore>1
            numTypeData(i) = 3;
        else
            numTypeData(i) = 1;
        end
    else
        if maxScore > 0
%             display(['altering ptftype interpretation: ' num2str(typeCounts(i,1)) ', ' num2str(typeCounts(i,2))]);
            numTypeData(i) = 0;
        else
            numTypeData(i) = 0;
        end
    end
    
    %ptfclass (e.g. supernova type ii/a)
    class = assetClassLabels(i);
    class = class{1};

    tokens = textscan(class, '%s');
    for t=1:length(tokens{1})
        token = tokens{1}{t};       
        if strcmp(token, '""')
            %do nothing wid dis            
        elseif ~isempty(strfind(token, 'SN')) || ...
                ~isempty(strfind(token, 'nova')) || ...
                ~isempty(strfind(token, 'SNI'))
            classCounts(i, 2) = classCounts(i, 2) + 1;  
        elseif ~isempty(strfind(token, 'AGN')) ||...%AGN = active galactic nucleus
                ~isempty(strfind(token, 'CV')) || ...%CV = cataclysmic variable star
                ~isempty(strfind(token, 'varstar')) ||...
                ~isempty(strfind(token, 'galaxy')) || ... 
                ~isempty(strfind(token, 'rock')) %might need to treat as unknown as rock classification is unreliable
            classCounts(i, 1) = classCounts(i, 1) + 1;
        elseif ~isempty(strfind(token, 'unknown'))
            classCounts(i, 3) = classCounts(i, 3) + 1;
        else
            %ignore other stuff but print it in case we missed something.
%             display(['Other string found in PTF class: ' token]);
        end
    end
    
    [maxScore, maxIdx] = max(classCounts(i,:));
    maxFlags = classCounts(i, :) == maxScore;
    if maxScore == 0
        numClassData(i) = 0;
    elseif sum(maxFlags) == 1
%         display(maxScore);
        if maxIdx==2
            numClassData(i) = 3;
        elseif maxIdx==1
            numClassData(i) = -1;
        elseif maxIdx==3
            numClassData(i) = 1;
        end
    else
        numClassData(i) = 1; 
    end
    
    if numTypeData(i) == 0 && maxScore > 0
        nClassOnly = nClassOnly + 1;
    end
end

%put the new values in place of the strings
assetTypeLabels = numTypeData;
assetClassLabels = numClassData;

knownTypeIdxs = assetTypeLabels~=0;
typeLabels = assetTypeLabels(knownTypeIdxs);
typeAssets = find(knownTypeIdxs);

knownClassIdxs = find(assetClassLabels~=0 & assetClassLabels~=1);
classLabels = assetClassLabels(knownClassIdxs);
% classAssets = assetIds(knownClassIdxs);

labels = sparse(knownClassIdxs, 1, double(classLabels), length(uniqueAssetIdxs), 1);
% labels = sparse(double(classAssets), ones(1, length(classAssets)), ...
%     double(classLabels), double(max(assets)), 1);
labels(labels==3) = 2;
labels(labels==-1) = 1;

%Treat "no type" as definitive no (used in papers)
if interpretNoTypes
    noTypes = assetTypeLabels==0 & assetClassLabels==0;
    labels(noTypes) = 1; 
end

%This is probably incorrect if the type labels are no wholly reliable -
%include these as reliable agents only
% %Currently, treat those with no from type labelling as definite no only
% negTypes = assetTypeLabels==-1 & assetClassLabels==0;
% labels(negTypes) = 1; 

snBaseOutputs = cell(3,1);
snBaseOutputs{1} = double(agentIdxs); %agents
snBaseOutputs{2} = double(assetIdxs); %assets
snBaseOutputs{3} = double(scores); %scores

typeAssets = double(typeAssets);
end