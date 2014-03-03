featureFile = '/homes/49/edwin/data/trec/finalOutput/MatrixTopic2000/Topic-Matrix/matrices_2000_topic.txt';
fileMapFile = '/homes/49/edwin/data/trec/finalOutput/MatrixTopic2000/Topic-Matrix/fileMap.txt';
labelFile = '/homes/49/edwin/data/trec/finalOutput/crowd_class_labels7.csv';

convertJudgementsToLabels

Tcols = dlmread(labelFile);

for l=1:size(Tcols,1)
    Tcols(l,3) = find(uniqueLabels==Tcols(l,3));
end

[workerIDs, ~, mappedWorkers] = unique(Tcols(:,1));

Tcols(:,1) = mappedWorkers;

Xworkers = {};
Xworkers{1} = Tcols(:,1);
Xworkers{2} = Tcols(:,2);
Xworkers{3} = Tcols(:,3);

Xfeat = dlmread(featureFile);

nDocs = Xfeat(1,1);
nFeat = Xfeat(2,1);
Xfeat = Xfeat(3:end,:);

threshold = 0;
Xfeat( Xfeat(:,3)<threshold, 3 ) = 0; 

%add one to indexes since matlab indexes from 1
Xfeat = sparse(Xfeat(:,2)+1, Xfeat(:,1)+1, Xfeat(:,3), nFeat, nDocs); 

nClasses = length(uniqueLabels);

Alpha0Work = zeros(nClasses, nClasses) + 0.1;
Alpha0Work(sub2ind([nClasses nClasses], 1:nClasses, 1:nClasses)) = 2;

Alpha0 = 2 .* ones(nClasses, 2); %0.5
Alpha0(:,2) = 1; % 0.5 0.25 is the best but very slow. Alternatives? Why is this so much better?


nClasses = size(Alpha0,1);

Nu0 = ones(1,nClasses);

bccSettings = settings.BccSettings();
featSettings = settings.BccSettings();
featSettings.Alpha = Alpha0;
bccSettings.nu{1} = Nu0;

bccSettings.changeRateMod = 1;% best 0.5
bccSettings.useLikelihood = false; %doesn't matter much
bccSettings.scoreMap = [];
bccSettings.IbccMapFunction = [];
bccSettings.debug = false;
bccSettings.convThreshold = 0.0001*nDocs;
bccSettings.convIt = 3;
bccSettings.maxIt = 100;
bccSettings.minScore = 1;
bccSettings.maxScore = 2;

nWorkers = length(workerIDs);

Xfeat = Xfeat ./ max(max(Xfeat,[],1));
Tvec = mappedLabelsSkipClashes;

%set an appropriate prior for the crowd's labels
Alpha0Feat = repmat(Alpha0, [1, 1, size(Xfeat,1)]);

bccSettings.Alpha = Alpha0Work;
featSettings.Alpha = Alpha0Feat;

combiner = combiners.bcc.MixedIbccVb(bccSettings, nFeat, featSettings, ...
    Xfeat, nFeat+nWorkers, 1, Tvec, [], nClasses, nClasses);
combiner.initAll = true;  
combinedPost = combiner.combineDecisions(Xworkers);


confMat = combiner.Alpha{1};
kappa = combiner.Nu ./ sum(combiner.Nu);

f1 = figure;
f2 = figure;

for w=1:nWorkers
    wIdxs = find(Tcols(:,1)==w);
    confw = confMat(:,:,wIdxs);
    
    acc = 0;
    for j=1:length(uniqueLabels)
        acc = acc +  kappa(j).* confw(j, j, :) ./ sum(confw(j,:,:), 2);%(kappa(j)./sum(kappa(2:end)))
    end
    acc = reshape(acc, 1, length(acc));
    set(0, 'CurrentFigure', f1);
    plot(acc);
    
    hold all
    
    set(0, 'CurrentFigure', f2);
    qt = combiner.q_t(wIdxs);
    plot(qt);
    
    hold all
end
