function [Pi_n, Alpha_n, filterMap] = GetPiSlice(Kt, n, nprev, Alpha)

Kfiltered = [];
NKfiltered = [];
if nargin < 4
    Kt = Kt(1:n);
    nprev = 0;
else
    Kt = Kt(nprev+1:n);
end
Kunique = unique(Kt); %the ids of the base classifiers

for k=Kunique'
    
    Nk = find(Kt==k); %the data points corresponding to this base classifier
    kExample = numel(Nk);%ceil(numel(Nk)/2); %take the middle of the slice so we get 
    NKfiltered = [NKfiltered Nk(kExample)];    
    Kfiltered = [Kfiltered k];
end
    
Alpha_n = Alpha(:,:,nprev+NKfiltered);

normTerm = sum(Alpha_n, 2);
normTerm = repmat(normTerm, 1, size(Alpha_n,2));
Pi_n = Alpha_n ./ normTerm;

filterMap = sparse(ones(length(NKfiltered),1), Kt(nprev+NKfiltered), (1:length(NKfiltered))'); %a map from the index in the complete set Kt