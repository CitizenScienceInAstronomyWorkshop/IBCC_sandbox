function Gmatrix = trust(labelFile, G)
% Assume turkers are labelled contiguously 1,2,3,....

data=load(labelFile);
if nargin < 2
    G=0.7; % Confidence G value (i.e. Pr(correct|"not confident")).
end

nts=max(data(:,1)); % Number of turkers.

validResps = find(data(:,5)~=0);
data = data(validResps,:);

cc=find(data(:,3)==data(:,5) & data(:,4)==1); % Correct and confident.
cn=find(data(:,3)==data(:,5) & data(:,4)==0); % Correct and not confident.
nc=find(data(:,3)~=data(:,5) & data(:,4)==1); % Incorrect and confident.
nn=find(data(:,3)~=data(:,5) & data(:,4)==0); % Incorrect and not confident.

% The n** are the numbers for each correctness/confidence pairs for
% each turker.
ncc=zeros(nts,1);
ncn=zeros(nts,1);
nnc=zeros(nts,1);
nnn=zeros(nts,1);

ncc = sparse( 1, data(cc,1), 1, 1, max(data(:,1)) );
ncn = sparse( 1, data(cn,1), 1, 1, max(data(:,1)) );
nnc = sparse( 1, data(nc,1), 1, 1, max(data(:,1)) );
nnn = sparse( 1, data(nn,1), 1, 1, max(data(:,1)) );

Gmatrix = (ncc+(1-G)*nnn+G*ncn+1)./(ncc+ncn+nnc+nnn+2);
Gmatrix = Gmatrix';
