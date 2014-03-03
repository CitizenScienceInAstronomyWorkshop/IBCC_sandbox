workerIG = dlmread('/homes/49/edwin/data/aaai/test_hiring/workerIG_partial.csv');

wIds = unique(workerIG(:,1));

x = zeros(length(wIds), size(workerIG,1));
ig = zeros(length(wIds), size(workerIG,1));
th = zeros(length(wIds), size(workerIG,1));
all = zeros(length(wIds), size(workerIG,1));

for l=1:size(workerIG,1)
    w = workerIG(l,1);
    k =find(wIds==w);
    i = find(x(k,:)==0,1);
    x(k, i) = l;
    ig(k, i) = workerIG(l,2);
    th(k,i) = workerIG(l,3);
    if l==1
        all(k,i) = workerIG(l,2)+workerIG(l+1,2)+workerIG(l+2,2) ./3;
    elseif l==size(workerIG,1)
        all(k,i) = workerIG(l-1,2)+workerIG(l,2)+workerIG(l-2,2) ./3;
    else
        all(k,i) = workerIG(l-1,2)+workerIG(l,2)+workerIG(l+1,2) ./3;
    end
end

% figure;
% for w=1:size(x,1)
%     stop = find(x(w,:)==0,1,'first') -1;
%     plot(x(w,1:stop), ig(w,1:stop));
%     hold all
% end

figure;

ig2 = zeros(length(wIds), 400);
th2 = zeros(length(wIds), 400);


for w=1:size(x,1)
    stop = find(x(w,:)==0,1,'first') -1;


    prev = 1;
    for i=1:stop
        ig2(w,prev:x(w,i)) = ig(w,i);
        th2(w,prev:x(w,i)) = th(w,i);

        prev = x(w,i)+1;
    end

end

meanIg = sum(ig2,1) ./ size(x,1);

for w=1:size(x,1)
    plot(ig2(w,:)-th2(w,:));
    hold all    
end