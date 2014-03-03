
commSizes = zeros(length(g),1);
for c=1:length(g)
    commSizes(c) = length(g{c});
end

[sorted sortIdxs] = sort(commSizes);

for c=sortIdxs'
    sizes = zeros(1, length(g{c}));
    for s=1:length(g{c})
        id = observedAgents(g{c}(s));
        idxs = find(testData{1}==id);
        sizes(s) = length(idxs);
    end
    
%     display([num2str(c) ': size=' num2str(length(g{c})) ', mean=' num2str(mean(sizes)) ', variance=' num2str(var(sizes))]); 
    display([num2str(c) ' & ' num2str(length(g{c})) ' & ' num2str(round(10*mean(sizes))/10) ' & ' num2str(round(10*var(sizes))/10) ' \\']); 
end