function SquaredHD = squaredHellinger(nAgents, nClasses, nScores, Pi, Kappa)
    %might need to add efficiency shortcuts in later.
    SquaredHD = ones(nAgents, nAgents);
    for j=1:nClasses
        for c=1:nScores
            PiVector = reshape(Pi(j, c, :), [1, nAgents]);
            SquaredHD = SquaredHD - (Kappa(j) .* (PiVector' * PiVector).^0.5);
        end
    end

    SquaredHD(SquaredHD<0) = 0;
end