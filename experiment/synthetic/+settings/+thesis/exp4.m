function [expSet, synthSet, bccSet, rerun, range, varName] = exp4(i)

    varName = 'Number of Inverted Agents';
    rerun = true; display('warning: likely to squash any saved datasets');
    settings.thesis.exp0;
    expSet.setDirNames([rootForAllExpts '/exp4_mvredo/']);

    %data generation
%     synthSet.pFlip_a = [0 0 0 0 0];
    range = 0:synthSet.nAgents;%[0 2 4 5 6 8 10];%:synthSet.nAgents;%[0 1 2 3 4 5];
%     if range(i) > 0
%         synthSet.pFlip_a(1:range(i)) = 1;
%     end
%     synthSet.means = [1 1 1 1 1; -1 -1 -1 -1 -1];
%     synthSet.deviations = [2 2 2 2 2; -2 -2 -2 -2 -2];
%     synthSet.cal = [1 1 1 1 1];

    synthSet.pMode = zeros(1, synthSet.nAgents);
    synthSet.pMode(1:range(i)) = 1;
end