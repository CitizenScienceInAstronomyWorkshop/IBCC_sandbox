function [ p_casej p_totalc p_total ] = targetCases2(k, X, nClasses, nScores, P)
%targetCases Calculate the probabilities of each combination of targets
%classified by k. 
%   Detailed explanation goes here

    X{2} = X{2}(X{1}==k);
    X{3} = X{3}(X{1}==k);
    
    nClassified = length(X{2});
        
    p_casej = zeros(nClasses, nScores, nClassified^2);
    p_totalc = zeros(nClasses, nScores, nClassified^2);
    p_total = zeros(nClasses, nScores, nClassified^2);
        
    
    %only accounts for static IBCC model, doesn't work out contribution of
    %old data points that has faded!!!!!
    %Can we detect influence by the amount of change in the smoothing step?
    %Not really. 
    %Is it also incorrect to treat the cases as independent?
    %Arse.
    for c=1:nScores
        for j=1:nClasses
            p_case = 1;
            p_totalj = 0;
            p_totaljc = 0;
            for n=1:nClassified %go through x? What about zeros????
                
                in = X{2}(n);
                
                p_case = reshape(p_case, numel(p_case), 1) * [P(in,j) 1-P(in,j)];
                p_totalj = repmat( reshape(p_totalj, numel(p_totalj), 1), 1, 2) + ...
                    repmat([1 0], numel(p_totalj),1);
                
                c_n = X{3}(n) == c;
                
                p_totaljc = repmat( reshape(p_totaljc, numel(p_totaljc), 1), 1, 2) + ...
                    repmat([c_n 0], numel(p_totaljc),1);
                
                %compress equivalent cases
                temp = p_totalj;
                
                p_case = sparse(temp+1, p_totaljc+1, p_case);
                [p_totalj, p_totaljc] = find(p_case);
                p_case = p_case(p_case>0);
                p_totalj = p_totalj-1;
                p_totaljc = p_totaljc-1;
                
%                 new_p_totalj = sparse(size(p_case, 1), size(p_case, 2) );
%                 idxs = sub2ind(size(new_p_totalj), temp+1, p_totaljc+1);
%                 new_p_totalj( idxs) =  p_totalj;
% 
%                 new_p_totaljc = sparse(size(p_case,1), size(p_case,2));
%                 new_p_totaljc( sub2ind(size(new_p_totaljc), temp+1, p_totaljc+1) ) =  p_totaljc;
% 
%                 p_totalj = new_p_totalj;
%                 p_totaljc = new_p_totaljc;
            end
            
            p_total(j,c,1:numel(p_totalj)) = reshape(p_totalj, [1,numel(p_totalj)]);
            p_totalc(j,c,1:numel(p_totalj)) = reshape(p_totaljc, [1,numel(p_totaljc)]);
            p_casej(j,c,1:numel(p_totalj)) = reshape(p_case, [1,numel(p_case)]);
        end
    end
end

