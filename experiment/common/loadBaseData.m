function [baseOutputs, rawData, labels, screeningLabels, screenedObjs, objIds] = ...
    loadBaseData(dataFilename, rawData, sepPosLabFile, sepNegLabFile, ...
        rawColumnMap, columnString, legalScores, ...
        interpretEmptyScreeningLabels, translateLabels)

%SETTINGS -----------------------------------------------------------------

%Outputs data from the raw file as the following:
%baseOutputs = a sparse representation (scoreset) of the base
%classifier outputs/features
%labels = the true labels (class labels in GZSN)
%screeningLabels = labels generated during a cheap screening process to weed
%out negative examples. The negative labels here should be accurate while the
%positive ones are not likely to be.
if nargin < 5 || isempty(rawColumnMap)
    %Map columns in the input CSV file to columns in baseOutputs and to
    %the label vectors.
    %The target indices are: 
    %1=classification ID
    %2=agent ID
    %3=assetID
    %4=score
    %5=screening label
    %6=confirmed label
    rawColumnMap = [1 2 3 6 4 5];
end

if nargin < 7 || isempty(legalScores)
    legalScores = [-1,1,3];
end

if nargin < 6 || isempty(columnString)
    columnString = '%u64 %u64 %u64 %s %s %d';
end

%Interpret missing screening labels as confirmed negative examples?
if nargin < 8
    interpretEmptyScreeningLabels = true;
end

if nargin < 9
    %default to the GZSN function as a template example
    translateLabels = @translateGZSNLabels;
end

% LOAD RAW DATA -----------------------------------------------------------

display('Loading data.');

if nargin > 2 && ~isempty(sepPosLabFile)
    fid = fopen(sepPosLabFile);
    sepPosLabs = textscan(fid, '%u64');
    sepPosLabs = sepPosLabs{1};
    fclose(fid);
else
    sepPosLabs = [];
end

if nargin > 3 && ~isempty(sepNegLabFile)
    fid = fopen(sepNegLabFile);
    sepNegLabs = textscan(fid, '%u64');
    sepNegLabs = sepNegLabs{1};
    fclose(fid);
else
    sepNegLabs = [];
end

if isempty(rawData)
    if strcmp(dataFilename(end-3:end), '.mat')
        load(dataFilename);
        if exist(snData,'var')
            %Supernovae data might be saved as variable snData
            rawData = snData;
        end
    else %assume .csv format
        display(['Reading input data from file ' dataFilename]);
        fid = fopen(dataFilename);
        rawData = textscan(fid, columnString, 'delimiter', ',');
        fclose(fid);
    end
end

%DUPLICATES ---------------------------------------------------------------

%Remove duplicates if we have the classification ID
cleanData = rawData;
if rawColumnMap(1)>0
    [cleanData{rawColumnMap(1)}, uniqueIdxs] = unique(rawData{rawColumnMap(1)});
    for i=1:length(cleanData)
        if i==rawColumnMap(1)
            continue %done it already
        end
        if ~isempty(rawData{i})
            cleanData{i} = rawData{i}(uniqueIdxs);
        end
    end
else
    cleanData = rawData;
end

agents = cleanData{rawColumnMap(2)};
objs = cleanData{rawColumnMap(3)};
scores = cleanData{rawColumnMap(4)};

%CLEAN IMPOSSIBLE SCORES --------------------------------------------------

%remove crappy data where score is impossible
validIdx = find(ismember(scores,legalScores));
scores = scores(validIdx);
objs = objs(validIdx);
agents = agents(validIdx);

[objIds uniqueObjIdxs objIdxs] = unique(objs);
[~, ~, agentIdxs] = unique(agents); 

if ~exist('sepPosLabs','var') || isempty(sepPosLabs)
    if rawColumnMap(5)>0
        screenFullVector = cleanData{rawColumnMap(5)};
        screenFullVector = screenFullVector(validIdx);
        screenFullVector = screenFullVector(uniqueObjIdxs);
    else
        screenFullVector = sparse(length(objIds),1);
    end
    
    if rawColumnMap(6)>0
        labelsFullVector = cleanData{rawColumnMap(6)};   
        labelsFullVector = labelsFullVector(validIdx);       
        labelsFullVector = labelsFullVector(uniqueObjIdxs);
    else
        labelsFullVector = sparse(length(objIds),1);
    end

    %APPLICATION-SPECIFIC LABEL HANDLING --------------------------------------

    %if the labels are different text fields that need to be converted to class
    %indexes, use an app-specific function here
    if ischar(labelsFullVector{1})
        [numericLabels, numericScreenLabels] = translateLabels(objIds, labelsFullVector, screenFullVector);
        %put the new values in place of the strings
        screenFullVector = numericScreenLabels;
        labelsFullVector = numericLabels;
    end
    
    labelVals = sort(unique(labelsFullVector));
    %0 is not a proper label it's a missing one
    if (length(labelVals)<3 && ismember(0,labelVals)) ...
         || (length(labelVals)<4 && ismember(0,labelVals) && ismember(-1,labelVals))
        labelsFullVector(labelsFullVector~=-1) = labelsFullVector(labelsFullVector~=-1) + 1;
    end
else
    %add the obj ids that were not in the crowdsourced list, if any
%     objIds = [objIds; sepPosLabs(~ismember(sepPosLabs,objIds))]; 
%     objIds = [objIds; sepNegLabs(~ismember(sepNegLabs,objIds))]; 
    
    labelsFullVector = zeros(1,length(objIds));
    labelsFullVector( ismember(objIds, sepPosLabs) ) = 2;
    labelsFullVector ( ismember(objIds, sepNegLabs) ) = 1;
    labelsFullVector ( ismember(objIds, sepPosLabs) & ismember(objIds, sepNegLabs) ) = 0;
       
    screenFullVector = sparse(1,length(objIds));
    interpretEmptyScreeningLabels = false;
end

% REFORMAT LABELS ---------------------------------------------------------

%screen labels translated to sparse representation
screenedIdxs = screenFullVector~=0;
screeningLabels = screenFullVector(screenedIdxs);
screenedObjs = find(screenedIdxs);

%class labels translated to sparse format where 0s were where there was no
%textual entry at all and -1 is where the entry suggests the class is unknown.
knownClassIdxs = find(labelsFullVector~=0 & labelsFullVector~=-1);
classLabels = labelsFullVector(knownClassIdxs);
% classAssets = objIds(knownClassIdxs);

%class labels translated to full length vector with sparse implementation
labels = sparse(knownClassIdxs, 1, double(classLabels), length(labelsFullVector), 1);
% labels = sparse(double(classAssets), ones(1, length(classAssets)), ...
%     double(classLabels), double(max(objs)), 1);

%Treat "no type" as definitive no (used in papers)
if interpretEmptyScreeningLabels
    noTypes = screenFullVector==0 & labelsFullVector==0;
    labels(noTypes) = 1; 
end
%This is probably incorrect if the type labels are not wholly reliable -
%include these as reliable agents only
negTypes = screenFullVector==-1 & labelsFullVector==0;
labels(negTypes) = 1; 

% FINALISE OUTPUTS --------------------------------------------------------

baseOutputs = cell(3,1);
baseOutputs{1} = double(agentIdxs); %agents
baseOutputs{2} = double(objIdxs); %objs
baseOutputs{3} = double(scores); %scores

%DUPLICATES AGAIN ---------------------------------------------------------
%remove duplicates manually
checkMat = sparse(baseOutputs{1}, baseOutputs{2}, 1);
checkMat( sub2ind(size(checkMat), baseOutputs{1}, baseOutputs{2}) ) = baseOutputs{3};
[I,J] = find(checkMat);

baseOutputs{1} = double(I);
baseOutputs{2} = double(J);
baseOutputs{3} = checkMat(sub2ind(size(checkMat),I,J));

screenedObjs = double(screenedObjs);
end