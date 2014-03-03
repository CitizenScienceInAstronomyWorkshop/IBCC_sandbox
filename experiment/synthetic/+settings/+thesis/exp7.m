
function [expSet, synthSet, bccSet, rerun, range, varName] = exp7(i)
    varName = 'Number of Known Labels';%: strength of over-exaggeration and under-exaggeration of confidence';
    rerun = false;
    settings.thesis.exp0;
    expSet.setDirNames([rootForAllExpts '/exp7_mvredo/']);
     
%     nStupid = 10;    
%     synthSet.nAgents = synthSet.nAgents + nStupid;

%     synthSet.cal = ones(1, synthSet.nAgents);
%     synthSet.cal(1:3) = [3 0.3 2];
%     synthSet.cal(6:10) = [2 3 0.3 0.5 0.7];
%     synthSet.bias = zeros(1, synthSet.nAgents);
    
%     synthSet.nInformedAgents = synthSet.nAgents;
%     synthSet.nInfSensors = synthSet.nInformedAgents;
%     synthSet.pMode = [zeros(1, synthSet.nAgents-nStupid) 0.5.*ones(1, nStupid)];    
        
    %add in 5 uninformative agents.
    nFlippers = 4;
    synthSet.nInformedAgents = synthSet.nAgents + nFlippers; %stupid agents are "informed" but with very error-prone sensors during testing, so get over-confident
    synthSet.nInfSensors = synthSet.nInfSensors + nFlippers;
    
    synthSet.pMode = [zeros(1, synthSet.nAgents) 0.5*ones(1, nFlippers)];    
    
    synthSet.nAgents = synthSet.nAgents + nFlippers;
    synthSet.cal = [synthSet.cal ones(1, nFlippers)];
    synthSet.bias = [synthSet.bias zeros(1, nFlippers)];
    synthSet.means = [synthSet.means repmat([1;-1],1,nFlippers)];
    synthSet.deviations = [synthSet.deviations repmat([1;-1],1,nFlippers)];
    
    %make sure we are not setting silly priors when most agents are stupid
%     bccSet.Alpha = [10 10; 10 10];
    
    expSet.propKnown = [0 0.002 0.005 0.01 0.015 0.025];
    range = expSet.propKnown.*expSet.nSamples;
    expSet.iPropKnown = i; 
end