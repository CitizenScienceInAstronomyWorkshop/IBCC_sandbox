figure; 

startAgent = 17;
nToPlot = 4;%size(Pi,3)

for a=startAgent:startAgent+nToPlot-1
    subplot(2, ceil(nToPlot/2), a-startAgent+1);
    X = 1:(size(Pi,1));
    bar3(X, Pi(:,:,a));    
    
%     set(gca, 'XTickLabel', 0:9);
%     set(gca, 'YTickLabel', 0:9);
%     
    xlabel('estimated class');
    ylabel('true class');
end