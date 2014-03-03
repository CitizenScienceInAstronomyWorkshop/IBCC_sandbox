%make averages of the Auc

aucs = aucs_dyn75;
x = x_75;

meanAucs = zeros(11, 800);

if iscell(aucs)
    for r=1:length(aucs)
        prevx = 0;
        for i=1:length(x{r})
            auc_i = aucs{r}(:, i);
            x_i = x{r}(i);
            if i==length(x{r})
                x_i = 800;
            end
            if x_i==0
                continue;
            end
            meanAucs(:, prevx+1:x_i)  = meanAucs(:,prevx+1:x_i) + repmat(auc_i, 1, x_i-prevx);
            prevx = x_i;
        end
    end
    meanAucs = meanAucs ./ length(aucs);
else
    prevx = 0;
    for i=1:length(x)
        auc_i = aucs(:, i);
        x_i = x(i);
        if x_i==0
            continue;
        end
        meanAucs(:, prevx+1:x_i)  = meanAucs(:,prevx+1:x_i) + repmat(auc_i, 1, x_i-prevx);
        prevx = x_i;
    end    
end

meanAucs = sum(meanAucs(2:end,:),1) ./ size(meanAucs(2:end,:),1);

stopIdx = find(meanAucs==0,1);
if ~isempty(stopIdx)
    meanAucs = meanAucs(1:stopIdx-1);
end
plot(meanAucs);