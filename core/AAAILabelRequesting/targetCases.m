function [ p_casej p_totalc p_total ] = targetCases(k, X, nClasses, nScores, P)
%targetCases Calculate the probabilities of each combination of targets
%classified by k. 
%   Detailed explanation goes here

    X{2} = X{2}(X{1}==k);
    X{3} = X{3}(X{1}==k);
    
    nClassified = length(X{2});
        
    p_casej = zeros(nClasses, nScores, 2^nClassified);
    p_totalc = zeros(nClasses, nScores, 2^nClassified);
    p_total = zeros(nClasses, nScores, 2^nClassified);
        
    for c=1:nScores
        
        for j=1:nClasses
            p_case = 1;
            p_totalj = 0;
            p_totaljc = 0;
            for n=1:nClassified %go through x?
                
                in = X{2}(n);
                
                p_case = reshape(p_case, numel(p_case), 1) * [P(in,j) 1-P(in,j)];
                p_totalj = repmat( reshape(p_totalj, numel(p_totalj), 1), 1, 2) + ...
                    repmat([1 0], numel(p_totalj),1);
                
                c_n = X{3}(n) == c;
                
                p_totaljc = repmat( reshape(p_totaljc, numel(p_totaljc), 1), 1, 2) + ...
                    repmat([c_n 0], numel(p_totaljc),1);
                
                %compress equivalent cases?                
            end
            
            p_total(j,c,:) = reshape(p_totalj, [1,1,numel(p_totalj)]);
            p_totalc(j,c,:) = reshape(p_totaljc, [1,1,numel(p_totaljc)]);
            p_casej(j,c,:) = reshape(p_case, [1,1,numel(p_case)]);
        end
    end


end

