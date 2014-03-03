function [labels] = generateLabels(nClasses, nSamples, expSettings)


    if nClasses==2
%         labels = binornd(ones(1,nSamples), expSettings.p_c1);
        labels = zeros(1,nSamples);
        labels(1:round(expSettings.p_c1*nSamples)) = 1;
    else
        display('Ignoring the p_c1 argument in expSettings. Just creating labels from uniform distribution');
        labels = unidrnd(nClasses, 1, nSamples) - 1;
    end
    
    dirName = datageneration.checkDataDir(expSettings);
    dlmwrite(sprintf('%s%s_labels.mat', dirName, expSettings.expLabel), labels);
end