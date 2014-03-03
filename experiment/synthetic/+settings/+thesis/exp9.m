function [expSet, synthSet, bccSet, rerun, range, varName] = exp9(i)

    varName = 'Sensor Error Rate (after agents are trained)';
    rerun = true; display('warning: likely to squash any saved datasets');
    settings.thesis.exp0;
    synthSet.nDatasets = 50;
    expSet.setDirNames('/homes/49/edwin/matlab/data/thesis_synth/exp9/');

    %data generation
    synthSet.means = [1 1 1 1 1; -1 -1 -1 -1 -1];

    %test with sensors all of same quality
    synthSet.deviations = [2 2 2 2 2; -2 -2 -2 -2 -2];

    synthSet.cal = [1 1 1 1 1];
    synthSet.bias = [0 0 0 0 0];

    range = [0 0.1 0.2 0.3 0.4 0.5];

    synthSet.pMode = range(i);
end