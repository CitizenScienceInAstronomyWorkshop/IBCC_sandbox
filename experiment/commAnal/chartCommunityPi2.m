function avgPis = chartCommunityPi2(g, P, Alpha, minCommSize)
%Bar charts for the Average Pi for a Community
figure; 
X = 1:(size(Alpha,1));
    
nCharts = 0;
if nargin < 4
    minCommSize = 0;
end

for c=1:numel(g)
    commAgents = g{c};
    
    if numel(commAgents) < minCommSize
        continue;
    else
        nCharts = nCharts + 1;
    end
    
end

nCharts

n = 0;
nRows = 4;

avgPis = zeros(size(Alpha,1), size(Alpha,2), length(g));

[membVals discreteComms] = max(P, [], 2);
for c=1:size(P,2)
    commAgents = find(discreteComms==c);
    if numel(commAgents) < minCommSize
        continue;
    else
        n = n + 1;
    end
        
    Alpha0 = [0.5 0.3 0.05; 0.18 0.36 0.41];
    AlphaComm = Alpha(:, :, commAgents) - repmat(Alpha0, [1, 1, length(commAgents)]);  %deduct the prior
    participation = reshape(P(commAgents, c), 1, 1, numel(commAgents));
    AlphaComm = AlphaComm .* repmat(participation, size(Alpha,1), size(Alpha,2));
    sumAlpha = sum(AlphaComm, 3) + Alpha0;
    normTerm = repmat(sum(sumAlpha, 2), 1, size(Alpha,2));
    avgPi = sumAlpha ./ normTerm;
    avgPis(:,:,c) = avgPi;

    subplot(nRows, ceil(nCharts/nRows), n);

    bar3(X, avgPi);    
    axis([0.5, 3.5, 0.5, 2.5, 0, 1]);
    set(gca,'XTickLabel',{'-1';'1';'3'})
    set(gca, 'YTickLabel', {'0'; '1'});
    set(gca,'FontSize', 10);
    if n==1
        xlabel('decision');
        ylabel('true class');
    end
    
    title(['Community ' num2str(c) ', ' num2str(numel(commAgents)) ' members']);
end
end