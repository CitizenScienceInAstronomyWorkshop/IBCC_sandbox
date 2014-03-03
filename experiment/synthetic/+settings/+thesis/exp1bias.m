function [expSet, synthSet, bccSet, rerun, range, varName] = exp1bias(i)

    varName = 'Consistent Bias of Agent Decisions';%: strength of over-exaggeration and under-exaggeration of confidence';
    rerun = true; display('warning: likely to squash any saved datasets');
    settings.thesis.exp0;
    expSet.setDirNames('/homes/49/edwin/matlab/data/thesis_synth/exp1/');

    %data generation
%     synthSet.means = [1 1 1 1 1; -1 -1 -1 -1 -1];

    %test with sensors all of same quality
%     synthSet.deviations = [1 1.5 3 9 10; 1 1.5 3 9 10];

%     synthSet.cal = [1 1 1 1 1];
    synthSet.bias = [0 0 0 0 0];

    range = [0 0.1 0.25 0.5 0.75];

    synthSet.bias = synthSet.bias + range(i);
end