M = length(aucs);
means = zeros(1,M);
stdDev = zeros(1,M);

for i=1:M
    means(i) = sum(aucs{i}) / length(aucs{i});
    stdDev(i) = std(aucs{i});
end