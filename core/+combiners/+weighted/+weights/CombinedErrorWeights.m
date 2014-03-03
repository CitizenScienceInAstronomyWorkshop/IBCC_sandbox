classdef CombinedErrorWeights < combiners.weighted.weights.Weights
    %PRECISIONWEIGHTS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        nAgents
        windowSize
        invW %W as calculated so far
        invW0 = 1; %initial weights
        beta = 0.5; %adjustment amount
        members
        nSamples
    end
    
    methods
        
        function obj = CombinedErrorWeights(nAgents, windowSize, members, nSamples)
            
            obj.nAgents = nAgents;
            obj.windowSize = windowSize;
            obj.nSamples = nSamples;
            
            obj.invW = sparse(nAgents, nSamples);
            %obj.W0 = ones(obj.nAgents, 1) ./ nAgents;
            obj.invW0 = nAgents .* ones(obj.nAgents, 1);
            obj.members = members;
        end
        
        function W = calculateWeights(obj, scores, post, T, combinedPost, n)               
            if n==1
                obj.invW = obj.invW0;
                W = 1./obj.invW;
                return;
            elseif n==0 || isempty(T)
                invWn = ones(obj.nAgents, 1) .* obj.nAgents;
                obj.invW = [obj.invW invWn];
                W = 1./obj.invW;
                return;
            end
            
            iStart = n - obj.windowSize;
            if iStart < 1
                  iStart = 1;
            end
            
            iEnd = n-1;
            if iEnd > size(T, 2)
                iEnd = size(T, 2);
            end
            
            %don't need to round - this is done by the combination function
            %if the combiner's output is binary. If not, we deal with
            %continous output of the combination function
            %rounding should be done already for maj voting
            %C = ones(obj.nAgents, 1) * round(combinedPost(iStart:n-1));
            
            %C = ones(obj.nAgents, 1) * combinedPost(iStart:iEnd);
            %T = ones(obj.nAgents, 1) * T(iStart:iEnd);
            T = T(:, iStart:iEnd);
            
            scores = scores(:, iStart:iEnd);
            scores(T==-1) = 0;
            combinedPost = combinedPost(:, iStart:iEnd);
            combinedPost(T(1,:)==-1) = 0;
            T(T==-1) = 0;
            
            %measure of agreement between combined and base learner
            %not using this at the moment - for continuous combined
            %outputs it may punish the wrong agents.
            %contribution = 1 - (abs(C - scores));
                        
            %replace contribution to mistakes with abs error:
            scores = sparse(abs(T - scores));
%             for i=1:iEnd-iStart
%                 scores(:,i) = abs(T(i) - scores(:,i));
%             end
            
            %don't count errors before the agent joined the cluster - we
            %will replace those with an estimate later
            if ~isempty(obj.members)
                
                M = obj.members(:, iStart:iEnd);
                %contribution = contribution .* M;
                scores = scores .*M;

                %if the cluster membership changes we have a problem. The overall
                %decision is not recalculated to take into account the effect
                %new agents would have made on past decisions. This makes
                %sense, as a change to the agent has probably occurred when it
                %changed membership.  New agents will have weight 1, which
                %would be higher than older members that have made a couple of
                %mistakes.
                %SOLUTION is to assume that the new agent is average for the
                %cluster; replace unknown contributions with average
                %contribution for samples with mistakes.

                %avg the errors of cluster members; 
                avgError = sum(scores, 1) ./ sum(obj.members(:, iStart:iEnd));

                %create a matrix of errors with nAgents rows.
                estimatedErrors = (ones(obj.nAgents, 1) * avgError);

                %apply the estimated errors for nonmembers 
                nonMembers = 1 - M;           
                estimatedErrors = estimatedErrors .* nonMembers;
                error = scores + estimatedErrors; %gives non members the average score
            else
                error = scores;
            end
            mistakes = error * abs(T(1,:) - combinedPost)';  
            %power = sum(mistakes, 2);
            %invWn = (obj.invW0 .^ power) ./ (obj.beta .^ mistakes);
            invWn = obj.invW0 ./ (obj.beta .^ mistakes);
            
            %set weights for agents outside the current cluster to zero
            if ~isempty(obj.members)
                invWn = invWn .* obj.members(:, n);
            end

            %normalise
            invWn = invWn ./ sum(invWn);           
            
            %expand the weight vector 
            %obj.invW = [obj.invW invWn];
            obj.invW(:, n) = invWn;
            invW = obj.invW;
            
            W = 1./ invW;
        end        
    end
    
end

