
figHandleInc = [];
figHandleDec = [];
figHandleIncdec = [];

for i=1:25
    load(['experiment_1_d' num2str(i) '_i1_combiners.mat']);
%     load(['experiment_2_d' num2str(i) '_i1_combiners.mat']);
    load('experiment_1_synthSet.mat');
%     load('experiment_2_synthSet.mat');
    
    Alpha = combs{3}.Alpha{1};
    Kt = combs{3}.C1_desparse;
    changePoints = synthSet.defNoiseChanges(i,:);
    
    [figHandleInc figHandleDec figHandleIncdec] = drawPiEvolution2Combined(Alpha, Kt, changePoints, figHandleInc, figHandleDec, figHandleIncdec, true);
end

%plot the static ones
figure

for i=1:25
    
    load(['experiment_1_d' num2str(i) '_i1_combiners.mat']);
    load('experiment_1_synthSet.mat');
    
    Alpha = combs{3}.Alpha{1};
    N = size(Alpha,3)./size(combs{2}.Alpha,3);
    
    Pi = combs{2}.Alpha ./ repmat(sum(combs{2}.Alpha, 2), 1, 2);

    for k=1:size(combs{2}.Alpha,3)
        patchline([0 N], [Pi(1,2,k) Pi(1,2,k)], 'EdgeColor', 'black', 'EdgeAlpha', 0.4);
        hold all
        patchline([0 N], [Pi(2,2,k) Pi(2,2,k)], 'EdgeColor', 'red', 'EdgeAlpha', 0.4);
        hold all
    end
end

axis([0 1 0 N]);
ylabel('Probability of Response=1');
xlabel('Sample no.');
title('E[\pi_j^k] for All Agents with Static IBCC');