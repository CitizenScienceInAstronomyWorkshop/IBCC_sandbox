function [expSet, synthSet, bccSet, rerun, range, varName] = exp3(i)

    varName = 'Number of Uninformative Agents';
    rerun = true; display('warning: likely to squash any saved datasets');
    settings.thesis.exp0;
    expSet.setDirNames([rootForAllExpts '/exp3_mvredo/']);

    %data generation
    range = [0 5 10 20 50];

    synthSet.means = [synthSet.means repmat([1;-1], 1, range(i))];
    synthSet.deviations = [synthSet.deviations repmat([1;-1], 1, range(i))];
    
    synthSet.nInformedAgents = synthSet.nAgents + range(i); %stupid agents are "informed" but with very error-prone sensors during testing, so get over-confident
    synthSet.nInfSensors = synthSet.nInfSensors + range(i);
    
    synthSet.pMode = [zeros(1, synthSet.nAgents) 0.5.*ones(1, range(i))];    
    
    synthSet.nAgents = synthSet.nAgents + range(i);
    synthSet.cal = [synthSet.cal ones(1, range(i))];
    synthSet.bias = [synthSet.bias zeros(1, range(i))];
    
    %make sure we are not setting silly priors when most agents are stupid
%     bccSet.Alpha = bccSet.Alpha + range(i) ./ (range(i)+max(max(bccSet.Alpha)) );
%     bccSet.Alpha = [10 10; 10 10];
end