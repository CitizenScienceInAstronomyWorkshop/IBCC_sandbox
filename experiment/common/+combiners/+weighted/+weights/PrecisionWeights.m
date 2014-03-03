classdef PrecisionWeights < combiners.weighted.weights.Weights
    %PRECISIONWEIGHTS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        agents
        newPrecision
        X %training data - needed if the agents do not supply precision themselves.
        
        members % record 0 and 1 for each agent at each sample to state whether we should ignore them or not
    end
    
    methods
        
        function obj = PrecisionWeights(agents, members, X, recalculate)
            obj.agents = agents; 
            obj.members = members;
            
            if isempty(agents) || isempty(agents{1}) || isempty(findprop(agents{1}, 'os'))...
                    || (nargin > 2 && recalculate)
                display('Will try to recalculate agent precisions.');
                obj.newPrecision = true;
                if ~exist('X','var') || isempty(X)
                    obj.X = [];
                else
                    obj.X = X;
                end
            else
                obj.newPrecision = false;
            end
        end
        
        function W = calculateWeights(obj, scores, post, T)

            nAgents = size(scores, 1); 
            
            if obj.newPrecision
                sigma = obj.calculateSD(scores, post, nAgents, T);
            else

                nSamples = length(scores);
                sigma = zeros(nAgents, nSamples);
                for a=1:nAgents
                    agent = obj.agents{a};
                    nSamplesTrained = size(agent.os.s, 1);
                    sigma(a, :) = agent.os.s(nSamplesTrained-nSamples+1:nSamplesTrained)';

                end
            end
            
            beta = 1.0 ./ (sigma .* sigma);
            
            if isempty(obj.members)
                obj.members = ones(size(post,1),size(post,2));
            end
            W = beta ./ (ones(nAgents, 1) * sum(beta .* obj.members, 1) );
        end
        
        function s = calculateSD(obj, logodds, post, nAgents, targets)
            
            nSamples = length(logodds);
            
            s = zeros(nAgents, nSamples);
            
            for a=1:nAgents
%                 if isempty(obj.agents)
%                     Xa = obj.X;
%                     dim = size(Xa, 2);
%                 else
%                     Xa = obj.X(:, find(obj.agents{a}.s));
%                     dim = obj.agents{a}.num_active_sensors;
%                 end
%                 
%                 I = eye(dim);
% 
%                 S = I;
%                 
%                 q=0;
                diffs = 0;
                for n=1:nSamples

                    %SD is the mean of the deviation from correct target
                    
                    
                    
                    if n > length(targets) || targets(n) == -1
                        diffs = diffs + post(a,n)*(1-post(a,n));
                    else
                        t = targets(n);
                        
                        diffs = diffs + (t - post(a,n))^2;
                    end
                    
                    s(a, n) = (diffs ./ n)  ^ 0.5;
                    
%                     x = Xa(n, :)';
%                     
%                     Sp = I*S*I' + q*I;
%                     s(a, n) = x' * Sp * x;
% 
%                     y = post(a, n);
%                     u = y*(1-y);
%                     kappa = (1 + (pi*(s(a, n) )/8)).^(-0.5);
%                     
%                     %calculate the covariance for this sample
%                     S = Sp - ((u * kappa^2) / (1 + u * s(a, n) *	kappa^2)) * (Sp * x) * ...
%                         (Sp * x)'; 
%                      
%                     %this is incorrect but not sure what else to do here
%                     if n < nSamples
%                         uup = post(a, n+1) * (1 - post(a, n+1));
%                     else
%                         uup = u;
%                     end
%                     q = max(uup - u, 0);
                end
            end
        end
             
    end
    
end

