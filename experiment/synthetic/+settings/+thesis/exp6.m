function [expSet, synthSet, bccSet, rerun, range, varName] = exp6(i)

    varName = 'Number of Agents that Duplicate Agent 1';
    rerun = true; display('warning: likely to squash any saved datasets');
    settings.thesis.exp0;
    expSet.setDirNames([rootForAllExpts '/exp6_mvredo/']);

    %data generation
%     synthSet.means = [1 1 1 1 1; -1 -1 -1 -1 -1];
% 
%     synthSet.deviations = synthSet.deviations .* 0.1;
%     synthSet.pMode = zeros(1,synthSet.nAgents) + 0.3;
    
    range = 0:(synthSet.nAgents+1);
    range = range;
    
    synthSet.cal = [synthSet.cal ones(1, range(i))];
    synthSet.bias = [synthSet.bias zeros(1, range(i))];
    synthSet.pMode = [0 0 0 0 0];% zeros(1, range(i))];
    synthSet.means = [synthSet.means repmat([1;-1], 1, range(i))];
    synthSet.deviations = [synthSet.deviations repmat([1;-1], 1, range(i))];
    
    synthSet.subs = agentexecution.subscribeSensors(synthSet.nAgents, synthSet.nSensors, synthSet);
    synthSet.trainingBatch = 1:synthSet.nAgents;
    for j=1:range(i)
        synthSet.subs{synthSet.nAgents+j} = synthSet.subs{1};
        synthSet.trainingBatch = [synthSet.trainingBatch 1];
    end
    
    synthSet.nAgents = synthSet.nAgents + range(i);
    synthSet.nInformedAgents =  synthSet.nAgents;
    
%     bccSet.Alpha = [1.1 1.0; 1.0 1.1];
end