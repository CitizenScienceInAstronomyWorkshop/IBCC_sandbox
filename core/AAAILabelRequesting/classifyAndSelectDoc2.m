function [samplesToRequest, workersToRequest, excludedWorkers, P, combiner, ...
    crowdLabels, origTsorted] = classifyAndSelectDoc2( nToTry, workerToAllocate, nToAllocate,...
    featureFile, labelFile, excludedWorkersFile, outputRequestFile, outputExcludeFile, ...
    topicDocNoPairsFile, crowdTrust, Alpha0, Nu0 )

%nRequests: number of new responses. Produce a new task for each worker
%that submitted a response, plus a task for unknown workers. This can cause
%some tasks to go stale waiting for a response, but allows us to more
%quickly and more optimally respond to active workers.
% A solution would be to re-assign all waiting tasks every n iterations

% display('classifying then selecting next document to label...');

%%%%%%%% Check Input Arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('Nu0','var') || isempty(Nu0)
    Nu0 = [1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000] .* 0.1;
end

if ~exist('nToTry','var') || isempty(nToTry)
    nToTry = 0;
end
if ischar(nToTry)
    nToTry = str2double(nToTry);
end

% if ~exist('nRequests','var')
%     nRequests = 1;
% elseif ischar(nRequests)
%     nRequests = str2double(nRequests);
% end
if ~exist('workerToAllocate','var')
    workerToAllocate = 0;
end

%%%%%%%%% Read feature file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xfeat = dlmread(featureFile); %should produce Xfeat

if size(Xfeat,2)==3 && max(max(Xfeat))>1 %sparse format - rejig this to be a matrix

    nDocs = Xfeat(1,1);
    nFeat = Xfeat(2,1);
    Xfeat = Xfeat(3:end,:);
        
    %add one to indexes since matlab indexes from 1
    Xfeat = sparse(Xfeat(:,1)+1, Xfeat(:,2)+1, Xfeat(:,3), nDocs,nFeat); 
else
    nFeat = size(Xfeat,2);
    nDocs = size(Xfeat,1);
end

%%%%%%%%%%%%%%% Preprocessing of input data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    crowdLabels = dlmread(labelFile);
catch exception
    display('Warning, exception while reading the label file. Probably just empty...');
    crowdLabels = zeros(0,5);
end

fid = fopen(topicDocNoPairsFile);
topicDocNoCells = textscan(fid, '%d %s', 'Delimiter', '," ', 'MultipleDelimsAsOne',true);

%list of possible classes
origTsorted = [0; sort(unique(topicDocNoCells{1}))]; 
nClasses = length(origTsorted);

if  ~isempty(crowdLabels) && (min(crowdLabels(:,3))>nClasses || min(crowdLabels(:,3))==0)
    %not yet converted crowdLabels to start from 1
    for d=1:length(crowdLabels(:,3))   
        %change topic IDs to indices starting from 1
        crowdLabels(d,3) = find(origTsorted==crowdLabels(d,3)); %assuming that the 0 values mean "none of the above"
        if size(crowdLabels, 2)>=5 && crowdLabels(d,5) > 0
            crowdLabels(d,5) = find(origTsorted==crowdLabels(d,5));
        end
    end
end

%%%%%%%%%%% Priors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('crowdTrust','var')
    crowdTrust = [10 11];
end

if ~exist('Alpha0','var') || isempty(Alpha0)
    Alpha0 = 0.5 .* ones(nClasses, 2); %0.5
    Alpha0(:,2) = 0.25; % 0.5 0.25 is the best but very slow. Alternatives? Why is this so much better?
end

bccSettings = settings.BccSettings();
featSettings = settings.BccSettings();
featSettings.Alpha = Alpha0;
bccSettings.nu{1} = Nu0;
bccSettings.changeRateMod = 1;% best 0.5
bccSettings.useLikelihood = false; %doesn't matter much
bccSettings.scoreMap = [];
bccSettings.IbccMapFunction = [];
bccSettings.debug = false;
bccSettings.convThreshold = 0.001*nDocs;
bccSettings.convIt = 3;
bccSettings.maxIt = 50;

%%%%%%%%%% Run Classifier %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output is P and posterior Alpha (should also output posterior Nu from IBCC)
    
workerMap = unique(crowdLabels(:,1));
for w=1:size(crowdLabels,1)
    crowdLabels(w,1) = find(crowdLabels(w,1)==workerMap);
end
nWorkers = length(workerMap);

for k=1:length(workerToAllocate)
    if ~ismember(workerToAllocate(k), workerMap)
        workerToAllocate(k) = 0;
    else
        workerToAllocate(k) = find(workerToAllocate(k)==workerMap);
    end
end

Xfeat = Xfeat ./ max(max(Xfeat,[],1));
Tvec = sparse(1, nDocs); %have no fully reliable labels except test items
if size(crowdLabels, 2)>=5
    Tvec(sub2ind(size(Tvec), ones(size(crowdLabels,1),1), crowdLabels(:,2))) = crowdLabels(:,5);
end

%turn the labels into an additional column for each class for each worker. 
% Ones in the column where the user has selected that class.
newCol = {}; 
newCol{1} = crowdLabels(:,1);
newCol{2} = crowdLabels(:,2);
newCol{3} = crowdLabels(:,3);

Xfeat = Xfeat';
Xworkers = newCol;

%set an appropriate prior for the crowd's labels
Alpha0Feat = repmat(Alpha0, [1, 1, size(Xfeat,1)]);
if size(crowdTrust,1)<nClasses
    Alpha0Work = zeros(nClasses, nClasses, 1) + crowdTrust(1);
    for j=1:nClasses
        Alpha0Work(j, j, :) = crowdTrust(2);
    end
else
    Alpha0Work = crowdTrust;
end
bccSettings.Alpha = Alpha0Work;
featSettings.Alpha = Alpha0Feat;

display('we should save and reload the classifier to get a head start each time');
display('could also reuse clusters on each iteration regardless of classifications?');
%cut sample size by removing useless docs
%possibly cut ig evaluation down to cadidates plus docs classified by agents

combiner = combiners.bcc.MixedIbccVb(bccSettings, nFeat, featSettings, ...
    Xfeat, nFeat+nWorkers, 1, Tvec, [], nClasses, nClasses);
combiner.initAll = true;            
[combinedPost Alpha] = combiner.combineDecisions(Xworkers);
P = combinedPost';

testResults = cell(size(P,2),1);
for j=1:nClasses
    testResults{j} = P(:, j)';
end

%%%%%%%%%% Label Selection Method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output is a list of samples (can now also produce a list of sample/worker pairs

if exist(excludedWorkersFile, 'file')
    try
        excludedWorkers = dlmread(excludedWorkersFile);
        for w=1:length(excludedWorkers)
            excludedWorkers(w) = find(excludedWorkers(w)==workerMap);
        end
    catch exception
%         display('Exception while reading excluded workers file - probably just empty');
        excludedWorkers = [];
    end
else
    excludedWorkers = [];
end

[excludedWorkers, samplesToRequest, workersToRequest] = ...
    selectIGPairs(workerToAllocate, nToTry, nToAllocate, Xworkers, nWorkers, P, ...
    combiner, excludedWorkers, bccSettings, nFeat, featSettings, Xfeat, Tvec, combiner.combinedPost, combiner.logJoint);

workerMap(nWorkers+1) = 0; %output a 0 for new workers
for t=1:length(workersToRequest)
    requestedIdx = workersToRequest(t);
    if requestedIdx~=0
        requestedIdx = workerMap(requestedIdx);
    end
    workersToRequest(t) = requestedIdx;
end
% display('Excluded workers: ');
% excludedWorkers
% display('WorkerMap');
% workerMap
% display('Including previously excluded workers...');
excludedWorkers = workerMap(excludedWorkers);

%Write task assignment requests to CSV file
fHandle = fopen(outputRequestFile, 'w+'); 
for t=1:length(workersToRequest)
    fprintf(fHandle, '%d,', samplesToRequest(t));
    fprintf(fHandle, '%d,', workersToRequest(t));
    fprintf(fHandle, '\n');
end
fclose(fHandle);

%Write excluded worker list to file
dlmwrite(outputExcludeFile, excludedWorkers);