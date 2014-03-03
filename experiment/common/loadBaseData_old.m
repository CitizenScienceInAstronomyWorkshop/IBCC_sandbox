function [baseOutputs, rawData, labels, screeningLabels, screenedObjs, objIds] = ...
    loadBaseData(dataFilename, rawData, rawColumnMap, legalScores, columnString, ...
        interpretEmptyScreeningLabels, translateLabels)

%SETTINGS -----------------------------------------------------------------

%Outputs data from the raw file as the following:
%baseOutputs = a sparse representation (scoreset) of the base
%classifier outputs/features
%labels = the true labels (class labels in GZSN)
%screeningLabels = labels generated during a cheap screening process to weed
%out negative examples. The negative labels here should be accurate while the
%positive ones are not likely to be.
if nargin < 3
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

if nargin < 4
    columnString = '%d %d %d %s %s %d';
end

if nargin < 5
    legalScores = [-1,1,3];
end

%Interpret missing screening labels as confirmed negative examples?
if nargin < 6
    interpretEmptyScreeningLabels = true;
end

if nargin < 7
    %default to the GZSN function as a template example
    translateLabels = @translateGZSNLabels;
end

% LOAD RAW DATA -----------------------------------------------------------

display('Loading data.');

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

%Remove duplicates
cleanData = rawData;
[cleanData{rawColumnMap(1)}, uniqueIdxs] = unique(rawData{rawColumnMap(1)});
for i=1:length(cleanData)
    if i==rawColumnMap(1)
        continue %done it already
    end
    if ~isempty(rawData{i})
        cleanData{i} = rawData{i}(uniqueIdxs);
    end
end
agents = cleanData{rawColumnMap(2)};
objs = cleanData{rawColumnMap(3)};
scores = cleanData{rawColumnMap(4)};
screenFullVector = cleanData{rawColumnMap(5)};
labelsFullVector = cleanData{rawColumnMap(6)};

%remove crappy data where score is impossible
validIdx = find(ismember(scores,legalScores));
scores = scores(validIdx);
objs = objs(validIdx);
agents = agents(validIdx);
screenFullVector = screenFullVector(validIdx);
labelsFullVector = labelsFullVector(validIdx);

%Ignore uncertain labels.
[objIds uniqueObjIdxs objIdxs] = unique(objs);
[agentIds uniqueAgentIdxs agentIdxs] = unique(agents); 
screenFullVector = screenFullVector(uniqueObjIdxs);
labelsFullVector = labelsFullVector(uniqueObjIdxs);

%APPLICATION-SPECIFIC LABEL HANDLING --------------------------------------

%if the labels are different text fields that need to be converted to class
%indexes, use an app-specific function here
if strcmp(class(screenFullVector{1}), 'char')
    [numericLabels, numericScreenLabels] = translateLabels(objIds, labelsFullVector, screenFullVector);
    %put the new values in place of the strings
    screenFullVector = numericScreenLabels;
    labelsFullVector = numericLabels;
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
labels = sparse(knownClassIdxs, 1, double(classLabels), length(uniqueObjIdxs), 1);
% labels = sparse(double(classAssets), ones(1, length(classAssets)), ...
%     double(classLabels), double(max(objs)), 1);

%Treat "no type" as definitive no (used in papers)
if interpretEmptyScreeningLabels
    noTypes = screenFullVector==0 & labelsFullVector==0;
    labels(noTypes) = 1; 
end
%This is probably incorrect if the type labels are not wholly reliable -
%include these as reliable agents only
% negTypes = screenFullVector==-1 & labelsFullVector==0;
% labels(negTypes) = 1; 

% FINALISE OUTPUTS --------------------------------------------------------

baseOutputs = cell(3,1);
baseOutputs{1} = double(agentIdxs); %agents
baseOutputs{2} = double(objIdxs); %objs
baseOutputs{3} = double(scores); %scores

screenedObjs = double(screenedObjs);
end