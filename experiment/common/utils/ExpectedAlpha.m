
%Alpha
%Need indices for agents we observed.

%use Alpha already calculated. 
%community detection

try
ELnPi(:, :, nonDefAgents) = psi(Alpha);
catch
    display('psi error');
end

normTerm = psi(sum(Alpha, 2));
normTerm = repmat(normTerm, 1, obj.nScores);
ELnPi(:,:,nonDefAgents) = ELnPi(:,:,nonDefAgents) - normTerm;