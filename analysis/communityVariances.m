for c=1:numel(gTasks1)
    
    mean = sum(Pi_relevant(:,:,gTasks1{c}), 3) ./ length(gTasks1{c});
    
    mean = repmat(mean, [1, 1, length(gTasks1{c})]);
    
    dists = sum((Pi_relevant(:,:,gTasks1{c})-mean).^2, 3) ./ (length(gTasks1{c}-1) );

    display(['community: ' num2str(c) ', ' num2str(length(gTasks1{c}))]);
    dists.^0.5
    
end