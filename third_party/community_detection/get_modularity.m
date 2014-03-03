function Q = get_modularity(partition,W,flag)
    if ~exist('flag','var')
        flag = 'groups';
    end

    if strcmp(flag,'groups')
        B = group_to_incidence_matrix(partition);
    else
        B = partition;
    end

    if size(B,2) == 1
        Q = 0;
    else
        e = (1/sum(sum(W)))* B' * W * B;

        Q = trace(e) - sum(sum(e^2));
    end

    function B = group_to_incidence_matrix(groups)
        K = length(groups);
        N = length(cat(2,groups{:}));

        B = zeros(N,K);

        for i=1:K
           B(groups{i},i) = 1; 
        end
    end
end