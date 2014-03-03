nClasses = 8;
nRequests = 200;

figure;

for j=1:nClasses
    subplot(2, nClasses/2, j);
    %limited to 300 because the later points are all basically the same
    data = aucs_combined(j,1:nRequests*5/3)';
    xpoints = [1:nRequests*5/3].*3;
    plot(xpoints, data);
    hold all
    
    data = aucs_knockon_5req(j,1:nRequests)';
    xpoints = [1:nRequests].*5;
    plot(xpoints, data);
    hold all

    data = aucs_mostuncertain(j,1:nRequests)'; 
    xpoints = [1:nRequests].*5;
    plot(xpoints, data);    
    hold all    
    
    data = aucs_random(j,1:nRequests)'; 
    xpoints = [1:nRequests].*5;
    plot(xpoints, data);    
    hold all
end