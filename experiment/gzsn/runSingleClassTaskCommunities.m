Alpha_static = agentRatings{1}{3};%runner.combiner{3}.Alpha;

class1Idxs = ismember(snBaseOutputs{2}, find(labels==1));
class2Idxs = ismember(snBaseOutputs{2}, find(labels==2));

class1Data = snBaseOutputs;
class1Data{1} = class1Data{1}(class1Idxs);
class1Data{2} = class1Data{2}(class1Idxs);
class1Data{3} = class1Data{3}(class1Idxs);

%static
[gTasks1_c1, Ptasks1_c1, gTasks2_c1, Ptasks2_c1, AlphaTaskComms_c1, PiTaskComms_c1, observedAgents_c1] ...
    = extractTaskCommunities(Alpha_static, class1Data, class1Data, labels);

[avgPis, gIdx] = chartCommunityPi(gTasks1_c1, Ptasks1_c1, AlphaTaskComms_c1, 1, 0);

validComms = find(sum(avgPis(1,:,:),2)>0);
avgPisCompact1 = avgPis(:,:,validComms);
gIdx1 = gIdx(validComms);
h = drawPiSlice(avgPisCompact1, gTasks1_c1, Ptasks1_c1, 1, gIdx1);

set(h, 'PaperPosition', [0 0 24 21.6]);
saveas(h, sprintf('/homes/49/edwin/Documents/BCC_paper_2012/new_graphs/commonTaskComm/no_subsampling/150312-2/%d_stat_c1.png', h), 'png');
saveas(h, sprintf('/homes/49/edwin/Documents/BCC_paper_2012/new_graphs/commonTaskComm/no_subsampling/150312-2/%d_stat_c1.fig', h), 'fig');
saveas(h, sprintf('/homes/49/edwin/Documents/BCC_paper_2012/new_graphs/commonTaskComm/no_subsampling/150312-2/%d_stat_c1.eps', h), 'epsc');

class2Data = snBaseOutputs;
class2Data{1} = class2Data{1}(class2Idxs);
class2Data{2} = class2Data{2}(class2Idxs);
class2Data{3} = class2Data{3}(class2Idxs);


%static
[gTasks1_c2, Ptasks1_c2, gTasks2_c2, Ptasks2_c2, AlphaTaskComms_c2, PiTaskComms_c2, observedAgents_c2] ...
    = extractTaskCommunities(Alpha_static, class2Data, class2Data, labels);

[avgPis, gIdx] = chartCommunityPi(gTasks1_c2, Ptasks1_c2, AlphaTaskComms_c2, 1, 0);

validComms = find(sum(avgPis(2,:,:),2)>0);
avgPisCompact2 = avgPis(:,:,validComms);
gIdx2 = gIdx(validComms);
h = drawPiSlice(avgPisCompact2, gTasks1_c2, Ptasks1_c2, 2, gIdx2);

set(h, 'PaperPosition', [0 0 24 21.6]);
saveas(h, sprintf('/homes/49/edwin/Documents/BCC_paper_2012/new_graphs/commonTaskComm/no_subsampling/150312-2/%d_stat_c2.png', h), 'png');
saveas(h, sprintf('/homes/49/edwin/Documents/BCC_paper_2012/new_graphs/commonTaskComm/no_subsampling/150312-2/%d_stat_c2.fig', h), 'fig');
saveas(h, sprintf('/homes/49/edwin/Documents/BCC_paper_2012/new_graphs/commonTaskComm/no_subsampling/150312-2/%d_stat_c2.eps', h), 'epsc');
