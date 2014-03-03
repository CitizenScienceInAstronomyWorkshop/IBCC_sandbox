function [avgPis, gIdx] = chartCommunityPi(g, P, Alpha, minCommSize, minObs, removePriors)
%Bar charts for the Average Pi for a Community
% figure; 
X = 1:(size(Alpha,1));

if nargin > 5 && removePriors==true
    Alpha0 = [0.5 0.3 0.05; 0.18 0.36 0.41];
    newAlpha0 = [0.05 0.05 0.05; 0.05 0.05 0.05]; %prevent us from getting silly ML results
    Alpha = Alpha(:, :, :) - repmat(Alpha0, [1, 1, size(Alpha,3)]) + repmat(newAlpha0, [1, 1, size(Alpha,3)]);  %deduct the prior    
end

if nargin < 5
    minObs = 0;
end

normTerm = sum(Alpha, 2);
normTerm = repmat(normTerm, 1, size(Alpha,2));
Pi = Alpha ./ normTerm;
    
nCharts = 0;
if nargin < 4
    minCommSize = 1;
end

for c=1:numel(g)
    commAgents = g{c};
    
    if numel(commAgents) < minCommSize
        continue;
    else
        nCharts = nCharts + 1;
    end
    
end

% nCharts

n = 0;
nRows = 4;

avgPis = zeros(size(Pi,1), size(Pi,2), length(g));

[membVals discreteComms] = max(P, [], 2);

gIdx = zeros(1, size(P,2));
gIdxCount = 1;

for c=1:size(P,2)
    commAgents = find(discreteComms==c);
%     commAgents = commAgents(P(commAgents, c)>0.5);
    if numel(commAgents) < minCommSize
        if numel(commAgents)>0
            gIdxCount = gIdxCount + 1;
        end
        continue;
    else
        n = n + 1;
        gIdx(c) = gIdxCount;
        gIdxCount = gIdxCount + 1;
    end
        
    PiComm = Pi(:, :, commAgents); 
    AlphaComm = Alpha(:, :, commAgents);
    
    PiTotals = repmat(sum(AlphaComm, 2), [1,size(Pi,2),1]);
    
    PiComm = PiComm .* repmat(reshape(P(commAgents, c), 1, 1, numel(commAgents)), size(Pi,1), size(Pi,2));
       
    PiComm(PiTotals<minObs) = 0;
   
    PiTotals(PiTotals<minObs) = 0;
    commSize = sum(PiTotals(:,1,:)~=0, 3);
%     for j=1:numel(commSize)
%         if commSize(j) == 0
%             PiTotals
%     end
    
    avgPi = sum(PiComm, 3);
        
%     avgPi = avgPi ./ sum(P(commAgents,c)); %sum(P(:,c), 1);
    %normalise
    avgPi = avgPi ./ repmat(sum(avgPi, 2), 1, size(Pi,2));
    avgPis(:,:,c) = avgPi;

%     subplot(nRows, ceil(nCharts/nRows), n);
% 
%     bar3(X, avgPi);    
%     axis([0.5, 3.5, 0.5, 2.5, 0, 1]);
%     set(gca,'XTickLabel',{'-1';'1';'3'})
%     set(gca, 'YTickLabel', {'0'; '1'});
%     set(gca,'FontSize', 10);
%     if n==1
%         xlabel('decision');
%         ylabel('true class');
%     end
end
% 
% title('Pi Means Weighted By Community Participation');
% 
% %draw prototypical individual only
% n=0;
% figure;
% [A, C] = max(P, [], 1);
% for c=1:size(P,2)
%     commAgents = find(discreteComms==c);
%     if numel(commAgents) < minCommSize
%         continue;
%     else
%         n = n + 1;
%     end
%     
%     agent = C(c)
%     
%     PiComm = Pi(:, :, agent)
%     subplot(nRows, ceil(nCharts/nRows), n);
% 
%     bar3(X, PiComm);     
% end
% 
% title('Prototypes (member with higest participation)');
end