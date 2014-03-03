% Alpha_static = agentRatings{1}{3};%runner.combiner{3}.Alpha;
% 
% %static
% [gTasks1, Ptasks1, gTasks2, Ptasks2, AlphaTaskComms, PiTaskComms, observedAgents] ...
%     = extractTaskCommunities(Alpha_static, testData, snBaseOutputs, labels);
% 
% % display('now with priors removed');
% [avgPis, gIdx] = chartCommunityPi(gTasks1, Ptasks1, AlphaTaskComms, 1, 1, true);
% 
% handles = [0 0];
% 
% validComms = find(sum(avgPis(1,:,:),2)>0);
% avgPisCompact1 = avgPis(:,:,validComms);
% gIdx1 = gIdx(validComms);
% handles(1) = drawPiSlice(avgPisCompact1, gTasks1, Ptasks1, 1, gIdx1);
% 
% validComms = find(sum(avgPis(2,:,:),2)>0);
% gIdx2 = gIdx(validComms);
% avgPisCompact2 = avgPis(:,:,validComms);
% handles(2) = drawPiSlice(avgPisCompact2, gTasks1, Ptasks1, 2, gIdx2);
% 
% for h=handles
%     set(h, 'PaperPosition', [0 0 24 21.6]);
%     saveas(h, sprintf('/homes/49/edwin/Documents/BCC_paper_2012/new_graphs/commonTaskComm/190312-4/%d_stat.png', h), 'png');
%     saveas(h, sprintf('/homes/49/edwin/Documents/BCC_paper_2012/new_graphs/commonTaskComm/190312-4/%d_stat.fig', h), 'fig');
%     saveas(h, sprintf('/homes/49/edwin/Documents/BCC_paper_2012/new_graphs/commonTaskComm/190312-4/%d_stat.eps', h), 'epsc');
% end
% close all
% % 
% % % display('now with min number of observations');
% % % avgPis = chartCommunityPi(gTasks1, Ptasks1, AlphaTaskComms, 1, 5);
% % % avgPisCompact = avgPis(:,:,find(sum(sum(avgPis))>0));
% % % drawPiSlice(avgPisCompact, gTasks1, Ptasks1, 1);
% % % drawPiSlice(avgPisCompact, gTasks1, Ptasks1, 2);
% % %  
% % %dynamic
% Alpha_all = agentRatings{1}{4}{1};%runner.combiner{4}.Alpha{1};
% % % % Pi_all = Alpha_all ./ repmat(sum(Alpha_all,2), 1, size(Alpha_all,2));%runner.combiner{4}.Pi;
% % % 
% slices = [1 100000 200000 493048];%[1 100000 200000 300000 400000 493048];%[1 25000 50000 100000 200000 300000 400000 493048];%[3400 5000 6640 13280 19920 26558];
% % slices = [1 3000 12000 26558];%[1 5300 10600 15900 21200 26558];[1 100000 300000 493048];
% gPiDyn = cell(1,numel(slices));
% PPiDyn = cell(1,numel(slices));
% gTasks1Dyn = cell(1,numel(slices));
% Ptasks1Dyn = cell(1,numel(slices));
% gTasks2Dyn = cell(1,numel(slices));
% Ptasks2Dyn = cell(1,numel(slices));
% 
% AlphaTaskCommsDyn = cell(1,numel(slices));
% PiTaskCommsDyn = cell(1,numel(slices));
% 
% Kt = testData{1};
% 
% handles1 = zeros(1, numel(slices)-1);
% handles2 = zeros(1, numel(slices)-1);
for s=2:numel(slices)
    n = slices(s);
    [Pi_slice, Alpha_slice, filterMap] = GetPiSlice(Kt, n, 0, Alpha_all);
%     [gTasks1Dyn{s-1}, Ptasks1Dyn{s-1}, gTasks2Dyn{s-1}, Ptasks2Dyn{s-1}, AlphaTaskCommsDyn{s-1}, PiTaskCommsDyn{s-1}, observedAgents] ...
%         = extractTaskCommunities(Alpha_slice, testData, snBaseOutputs, labels, slices(s), filterMap);%, slices(s-1));    
%     [gPiDyn{s}, PPiDyn{s}] ...
%         = extractCommunities(Alpha_slice, testData, snBaseOutputs, true);
end
% 
% 
for s=1:numel(slices)-1
% %     display(['now dynamic for slice ' num2str(slices(s))]);
% %     [avgPis, gIdx] = chartCommunityPi(gTasks1Dyn{s}, Ptasks1Dyn{s}, AlphaTaskCommsDyn{s}, 1, 2);
% %     avgPisCompact = avgPis(:,:,find(sum(sum(avgPis))>0));
% %     gIdx = gIdx(find(sum(sum(avgPis))>0));    
% %     drawPiSlice(avgPisCompact, gTasks1Dyn{s}, Ptasks1Dyn{s}, 1, gIdx);
% %     drawPiSlice(avgPisCompact, gTasks1Dyn{s}, Ptasks1Dyn{s}, 2, gIdx);
%     
%     display('now with priors removed');
    [avgPis, gIdx] = chartCommunityPi(gTasks1Dyn{s}, Ptasks1Dyn{s}, AlphaTaskCommsDyn{s}, 1, 1, false);
%     
%     validComms = find(sum(avgPis(1,:,:),2)>0);
%     gIdx1 = gIdx(validComms);
%     avgPisCompact1 = avgPis(:,:,validComms);    
%     handles1(s) = drawPiSlice(avgPisCompact1, gTasks1Dyn{s}, Ptasks1Dyn{s}, 1, gIdx1);
%     
%     validComms = find(sum(avgPis(2,:,:),2)>0);
%     gIdx2 = gIdx(validComms);   
%     avgPisCompact2 = avgPis(:,:,validComms);
%     handles2(s) = drawPiSlice(avgPisCompact2, gTasks1Dyn{s}, Ptasks1Dyn{s}, 2, gIdx2);
%     
% %     display('now with min number of observations'); %don't use: not
% %     enough observations with the subsampled dataset
% %     avgPis = chartCommunityPi(gTasks1Dyn{s}, Ptasks1Dyn{s}, AlphaTaskCommsDyn{s}, 1, 5);
% %     avgPisCompact = avgPis(:,:,find(sum(sum(avgPis))>0));
% %     drawPiSlice(avgPisCompact, gTasks1Dyn{s}, Ptasks1Dyn{s}, 1);
% %     drawPiSlice(avgPisCompact, gTasks1Dyn{s}, Ptasks1Dyn{s}, 2);
end
% 
% for h=handles1
%     set(h, 'PaperPosition', [0 0 24 21.6]);
%     saveas(h, sprintf('/homes/49/edwin/Documents/BCC_paper_2012/new_graphs/commonTaskComm/190312-4/alpha1_%d_dyn1.png', h), 'png');
%     saveas(h, sprintf('/homes/49/edwin/Documents/BCC_paper_2012/new_graphs/commonTaskComm/190312-4/alpha1_%d_dyn1.fig', h), 'fig');
%     saveas(h, sprintf('/homes/49/edwin/Documents/BCC_paper_2012/new_graphs/commonTaskComm/190312-4/alpha1_%d_dyn1.eps', h), 'epsc');
% end
% 
% for h=handles2
%     set(h, 'PaperPosition', [0 0 24 21.6]);
%     saveas(h, sprintf('/homes/49/edwin/Documents/BCC_paper_2012/new_graphs/commonTaskComm/190312-4/alpha1_%d_dyn2.png', h), 'png');
%     saveas(h, sprintf('/homes/49/edwin/Documents/BCC_paper_2012/new_graphs/commonTaskComm/190312-4/alpha1_%d_dyn2.fig', h), 'fig');
%     saveas(h, sprintf('/homes/49/edwin/Documents/BCC_paper_2012/new_graphs/commonTaskComm/190312-4/alpha1_%d_dyn2.eps', h), 'epsc');
% end

% display('press a key to continue drawing the common task-score means');
% pause;
% 
% avgPis = chartCommunityPi(gTasks2, Ptasks2, AlphaTaskComms, 1);
% avgPisCompact = avgPis(:,:,find(sum(sum(avgPis))>0));
% drawPiSlice(avgPisCompact, gTasks2, Ptasks2, 1);
% drawPiSlice(avgPisCompact, gTasks2, Ptasks2, 2);
% 
% % display('now with priors removed');
% % avgPis = chartCommunityPi(gTasks2, Ptasks2, AlphaTaskComms, 1, 1, true);
% % avgPisCompact = avgPis(:,:,find(sum(sum(avgPis))>0));
% % drawPiSlice(avgPisCompact, gTasks2, Ptasks2, 1);
% % drawPiSlice(avgPisCompact, gTasks2, Ptasks2, 2);
% 
% % display('now with min number of observations');
% % avgPis = chartCommunityPi(gTasks2, Ptasks2, AlphaTaskComms, 1, 5);
% % avgPisCompact = avgPis(:,:,find(sum(sum(avgPis))>0));
% % drawPiSlice(avgPisCompact, gTasks2, Ptasks2, 1);
% % drawPiSlice(avgPisCompact, gTasks2, Ptasks2, 2);
% 
% for s=1:numel(slices)
%     display(['now dynamic for slice ' num2str(slices(s))]);    
%     avgPis = chartCommunityPi(gTasks2Dyn{s}, Ptasks2Dyn{s}, AlphaTaskCommsDyn{s}, 1);
%     avgPisCompact = avgPis(:,:,find(sum(sum(avgPis))>0));
%     drawPiSlice(avgPisCompact, gTasks2Dyn{s}, Ptasks2Dyn{s}, 1);
%     drawPiSlice(avgPisCompact, gTasks2Dyn{s}, Ptasks2Dyn{s}, 2);
%     
% %     display('now with priors removed');
% %     avgPis = chartCommunityPi(gTasks2Dyn{s}, Ptasks2Dyn{s}, AlphaTaskCommsDyn{s}, 1, 1, true);
% %     avgPisCompact = avgPis(:,:,find(sum(sum(avgPis))>0));
% %     drawPiSlice(avgPisCompact, gTasks2Dyn{s}, Ptasks2Dyn{s}, 1);
% %     drawPiSlice(avgPisCompact, gTasks2Dyn{s}, Ptasks2Dyn{s}, 2);
%     
% %     display('now with min number of observations');
% %     avgPis = chartCommunityPi(gTasks2Dyn{s}, Ptasks2Dyn{s}, AlphaTaskCommsDyn{s}, 1, 5);
% %     avgPisCompact = avgPis(:,:,find(sum(sum(avgPis))>0));
% %     drawPiSlice(avgPisCompact, gTasks2Dyn{s}, Ptasks2Dyn{s}, 1);
% %     drawPiSlice(avgPisCompact, gTasks2Dyn{s}, Ptasks2Dyn{s}, 2);
% end
