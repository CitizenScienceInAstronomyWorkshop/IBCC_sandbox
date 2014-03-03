function [expSet, synthSet, bccSet, rerun, range, varName] = exp8(i)

    varName = 'Sensor Error for Agents 1-3';

%     varName = 'Miscalibration of Agents';%: strength of over-exaggeration and under-exaggeration of confidence';
    rerun = true; display('warning: likely to squash any saved datasets');
    settings.thesis.exp0;
    expSet.setDirNames('/homes/49/edwin/matlab/data/thesis_synth/exp8/');

    %data generation
%     synthSet.means = [1 1 1 1 1; -1 -1 -1 -1 -1];

    %test with sensors all of same quality
%     synthSet.deviations = [1 1.5 3 9 10; 1 1.5 3 9 10];
%     synthSet.deviations = [2 2 2 2 2; -2 -2 -2 -2 -2];

%     synthSet.cal = [5 5 2.5 0.5 0.5];
%     synthSet.cal = [2 2 0.5 0.5 1];
%     synthSet.bias = [0 0 0 0 0];
    
    range = [0 0.25 0.5 1 2 2.5 3];

    synthSet.cal = synthSet.cal .^ range(i);
end