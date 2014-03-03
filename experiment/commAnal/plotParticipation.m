figure;
m=1;
n=4;

K = [79142 142372 139963 259297];
for i=1:length(K)
    subplot(m,n,i)
    
    k = K(i);
    
    data = zeros(3, 5); 
    sliceId = filterMapSlices{3}(k);
    vals = PSlices{3}(sliceId,:);
    data(1, 1:length(vals)) = vals;
    
    sliceId = filterMapSlices{3}(k);
    vals = PSlices{12}(sliceId,:);
    data(2, 1:length(vals)) = vals;
    
    sliceId = filterMapSlices{3}(k);
    vals = PSlices{27}(sliceId,:);
    data(3, 1:length(vals)) = vals;
    
    bar3(data);
    
    xlabel('community');
    set(gca, 'YTickLabel', [3000 12000 26558]);
    ylabel('no. observations');
    title(['Base Classifier ' num2str(k)]);
end