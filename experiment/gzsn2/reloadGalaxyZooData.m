function [snBaseOutputs, labels, typeLabels, classLabels, typeOnlyLabels, ...
    typeAssets, classAssets] = reloadGZSNData(loadDataFromFile, dataIsAnnotations)
global snRawData 

if loadDataFromFile && isempty(snRawData)
    %%load('/home/edwin/work_overflow/galaxyZoo3/snData.mat');
    load('/homes/49/edwin/matlab/combination/data/galaxyZoo3/snData.mat');
    snRawData = snData;
else
    snData = snRawData;
end

%If we are using a full list of annotations use this version
if nargin > 1 && dataIsAnnotations
    %classifications only, not the individual annotations
    classData = snData;

    [classData{2}, uniqueIdxs] = unique(snData{2});
    for i=1:length(classData)
        if i==2
            continue
        end
        classData{i} = snData{i}(uniqueIdxs);
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
validIdx = find(scores~=2 & scores~=0);
scores = scores(validIdx);
assets = assets(validIdx);
agents = agents(validIdx);
ptfType = ptfType(validIdx);
ptfClass = ptfClass(validIdx);

assetIds = assets;
[vals idxs] = unique(assetIds);
assetIds = assets(idxs); %get the unique asset IDs - previously they were duplicated across each zoo classification
assetTypeLabels = ptfType(idxs);
assetClassLabels = ptfClass(idxs);

%Translate the PTF labels into scores. 
%At the moment we have all base classifiers emitting scores of the same
%type and from the same set. This allows only -1, 1 and 3.
%ptftype (e.g. transient)
% map: non-transient -> -1
% transient -> 1
% more than one transient labelling -> 3
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
%         tokenStr = token{1};
        if ~isempty(strfind(token, 'Transient'))
            typeCounts(i, 2) = typeCounts(i, 2) + 1;
        else
            typeCounts(i, 1) = typeCounts(i, 1) + 1;
        end
    end
    [maxScore, maxIdx] = max(typeCounts(i,:));
    maxFlags = typeCounts(i, :) == maxScore;
    if sum(maxFlags) == 1
        if maxIdx==1
            numTypeData(i) = -1;
        elseif maxScore>1
            numTypeData(i) = 3;
        else
            numTypeData(i) = 1;
        end
    else
        numTypeData(i) = 0;
    end
    
%ptfclass (e.g. supernova type ii/a)
    class = assetClassLabels(i);
    class = class{1};
    %AGN = active galactic nucleus
    %CV = cataclysmic variable star
    tokens = textscan(class, '%s');
    for t=1:length(tokens{1})
        token = tokens{1}{t};       
        if strcmp(token, '""')
            %do nothing wid dis            
        elseif ~isempty(strfind(token, 'SN')) || ...
                ~isempty(strfind(token, 'nova')) || ...
                ~isempty(strfind(token, 'SNI'))
            classCounts(i, 2) = classCounts(i, 2) + 1;  
        elseif ~isempty(strfind(token, 'AGN')) ||...
                ~isempty(strfind(token, 'CV')) || ....
                ~isempty(strfind(token, 'varstar')) ||...
                ~isempty(strfind(token, 'galaxy'))
            classCounts(i, 1) = classCounts(i, 1) + 1;
        elseif ~isempty(strfind(token, 'unknown')) ||...
                ~isempty(strfind(token, 'rock'))
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
        %display(maxScore);
        if maxIdx==2
%             display([type ', ' class]);
            numClassData(i) = 3;
        elseif maxIdx==1
            numClassData(i) = -1;
        elseif maxIdx==3
            numClassData(i) = 1;
        end
    else
        numClassData(i) = 1; 
        %display(maxFlags);
    end
    
    if numTypeData(i) == 0 && maxScore > 0
        nClassOnly = nClassOnly + 1;
    end
    
    %display(i);
end

%put the new values in place of the strings
assetTypeLabels = numTypeData;
assetClassLabels = numClassData;

% assetClassLabels = assetTypeLabels; %see if classification into types is any better

knownTypeIdxs = assetTypeLabels~=0;
typeLabels = assetTypeLabels(knownTypeIdxs);
typeAssets = assetIds(knownTypeIdxs);

knownClassIdxs = assetClassLabels~=0 & assetClassLabels~=1;
% knownClassIdxs = assetClassLabels ~= 0;
classLabels = assetClassLabels(knownClassIdxs);
classAssets = assetIds(knownClassIdxs);

labels = sparse(double(classAssets), ones(1, length(classAssets)), ...
    double(classLabels), double(max(assets)), 1);
% labels(labels==1) = 0;%1.5; %we will fix this earlier so that classLabels
% and classAssets do not contain unknown labels
labels(labels==3) = 2;
labels(labels==-1) = 1;

% labels(labels==1) = 2;
% labels(labels==-1) = 1;

knownTypeOnlyIdxs = knownTypeIdxs - knownClassIdxs;
knownTypeOnlyIdxs = knownTypeOnlyIdxs>0;
typeOnlyLabels = assetTypeLabels(knownTypeOnlyIdxs);
typeOnlyAssets = assetIds(knownTypeOnlyIdxs);

%These are not definitive so don't use as binary labels
% negTypes = find(typeLabels==-1);
% negTypeAssets = typeAssets(negTypes);
% updatableLabels = find(labels(negTypeAssets)==0);
% labels(negTypeAssets(updatableLabels)) = 1;

%Treat "no type" as definitive no
noTypes = assetIds(assetTypeLabels==0 & assetClassLabels==0);
labels(noTypes) = 1; 
%a bunch more assets have a type or class label that is no definite so we don't use it as a label

assetClassText = ptfClass(idxs); %text for each unique asset
classText = assetClassText(knownClassIdxs); %text for each asset that is given a class label of SN or not SN


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
snBaseOutputs = cell(3,1);
snBaseOutputs{1} = double(agents); %agents
snBaseOutputs{2} = double(assets); %assets
snBaseOutputs{3} = double(scores); %scores



end