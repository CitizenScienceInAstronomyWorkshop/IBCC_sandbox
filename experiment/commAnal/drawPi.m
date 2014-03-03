function drawPi(Alpha)

normTerm = sum(Alpha, 2);
normTerm = repmat(normTerm, 1, size(Alpha,2));
Pi = Alpha ./ normTerm;

figure;

nRows = 4;

for a=1:size(Pi,3)

    PiComm = Pi(:, :, a);
    subplot(nRows, ceil(size(Pi,3)/nRows), a);
    
    bar3(PiComm);
end
end