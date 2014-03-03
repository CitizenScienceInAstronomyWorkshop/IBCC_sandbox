function [Pi] = recoverPi(agentRatings, foldIdx, methodIdx)

Alpha = agentRatings{foldIdx}{methodIdx};

nScores = size(Alpha,2);
 
normTerm = sum(Alpha, 2);
normTerm = repmat(normTerm, 1, nScores);
Pi = Alpha ./ normTerm;
end