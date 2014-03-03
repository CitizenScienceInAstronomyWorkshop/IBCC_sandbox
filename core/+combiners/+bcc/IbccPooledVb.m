classdef  IbccPooledVb < combiners.bcc.IbccVb
 
    properties 
	%blending parameter that controls the mixing between the data for a row in a confusion matrix 
	%and the pooled vector; it also controls the way that the pooled vector 
	%weights data from each row in the matrix
	B = 20;%15;%value used in the emailed results %referred to as alpha in Tom Leonard 1977 renamed to avoid confusion.
	
	beta;
	gamma;
    end
 
    methods (Static)
        function sl = shortLabel
            sl = 'IbccPooledVb';
        end   
    end
    
    methods       
        function obj = IbccPooledVb(bccSettings, nAgents, K, targets, agents, nClasses, nScores)
            obj@combiners.bcc.IbccVb(bccSettings, nAgents, K, targets, agents, nClasses, nScores);
            obj.label = 'VB I.B.C.C. with pooling';
            
            obj.beta = obj.nScores*0.1;
            obj.gamma = repmat([0.6 0.3 0.1], [1, 1, obj.nAgents]);
            obj.settings.priorAdjLevel = 0; %ignore this bit of fiddling - this class is intended to avoid it completely
        end
                
        function initAlpha(obj, Alpha)
            if isempty(obj.Alpha0)
                obj.Alpha0 = zeros(obj.nClasses, obj.nScores, obj.nAgents);

                obj.Alpha0(:, :, :) = 1 ./ obj.nScores;      
                
                %sequential version needs stronger priors to avoid assuming too much from the early data points
                if obj.sequential
                    obj.Alpha0 = Alpha + 2;
                end
            else
                if size(obj.Alpha0, 3) < obj.nAgents
                    obj.setAlphaPrior(obj.Alpha0, true, obj.nAgents);
                end
            end       
%             for j=1:obj.nClasses
%                 level = obj.settings.priorAdjLevel;
%                 agentCounts = full(sum(obj.Cmat,2));
% 
%                 changeIdxs = agentCounts<level & agentCounts>0;
%                 
%                 priorScale = reshape(repmat( level./agentCounts(changeIdxs), obj.nScores, 1), ...
%                     [1, obj.nScores, sum(changeIdxs)]);
%                 obj.Alpha0(j,:,changeIdxs) = obj.Alpha0(j,:,changeIdxs) .* priorScale;
%             end
        end
        
        function Xi = updateApproxXi(obj, Count, ET)
	    %rowTotals = repmat(sum(Count,3), [1, 1, obj.nAgents]);
	    rowTotals = repmat(sum(Count,2), 1, obj.nScores);
	    tau = (1+obj.B) ./ (rowTotals + obj.B);
	    
	    Xi = (sum(tau.*Count,1) + obj.beta.*obj.gamma) ./ (sum(tau.*rowTotals,1) + obj.beta);
	    Xi = repmat(Xi, obj.nClasses, 1);	    
        end
                
        function Alpha = updateAlpha(obj, C, ET, agentIdx, OldPostAlpha, OldPostAlpha_t)
            if nargin < 5
                Count = obj.voteCounts(C, ET);
                Alpha = obj.Alpha0 + Count;
                Alpha = Alpha + obj.B .* obj.updateApproxXi(Alpha,ET);
            else
                Count_a = obj.voteCounts(C, ET, agentIdx);
                Alpha = OldPostAlpha; %modified from tested version which appeared to be wrong
                Alpha(:,:,agentIdx) = obj.Alpha0(:,:,agentIdx) + Count_a;
            end                 
        end
    end
end