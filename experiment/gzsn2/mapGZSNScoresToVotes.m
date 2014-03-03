function [ mappedC ] = mapGZSNScoresToVotes( C )
%     C(C>1) = 1;
%     mappedC = round(C/2 + 1.5);
    mappedC = (C+3)/2;
end

