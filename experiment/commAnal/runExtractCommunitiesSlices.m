PiSlices = cell(1,27);
AlphaSlices = cell(1,27);
filterMapSlices = cell(1,27);

selectedSlices = [6 10 14 18 22];

for i=selectedSlices
    [PiSlices{i} AlphaSlices{i} filterMapSlices{i}] = GetPiSlice(testData{1}, i*1000, 0, Alpha);
end

[PiSlices{27} AlphaSlices{27} filterMapSlices{27}] = GetPiSlice(testData{1}, 26558, 0, Alpha);

gSlices = cell(1,27);
PSlices = cell(1,27);
avgPiSlices = cell(1,27);

for i=selectedSlices
    reply = 'y';
    
    while strcmp('y', reply)
        [gSlices{i} PSlices{i} avgPiSlices{i}] = extractCommunities(AlphaSlices{i}, testData, snBaseOutputs, true);
        display([num2str(i) ', ' num2str(numel(gSlices{i}))]);
        reply = input('repeat?' , 's');
    end    
end