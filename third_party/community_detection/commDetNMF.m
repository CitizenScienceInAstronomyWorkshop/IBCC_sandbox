function [P,g,W,H] = commDetNMF(V,max_rank,W0,H0)
% Implements the community detection methodology presented in 
% "Overlapping Community Detection using Nonnegative Matrix Factorization", Physical Review E 83, 066114 (2011)
% by Ioannis Psorakis, Stephen Roberts, Mark Ebden and Ben Sheldon.
% University of Oxford, United Kingdom.

% Inputs :
% V: NxN adjacency matrix
% max_rank: a prior estimation on the maximum number of communities
% (optional).

% Outputs :
% 
% P: NxK matrix where each element Pik represents the (normalised)
% soft membership score of node-i to each community-k. P captures the OVERLAPPING
% community structure.

% g: groups cell array. Each element g{k} contains the node indices of
% community k. Each node i is placed in the community g{k} for which it has
% the maximum Pnk. Therefore g describes the NON-OVERLAPPING community
% structure.
%
% W,H are the NMF factors (normalised from 0 to 1)

% NOTES:
% ------
% - In many cases the algorithm gives better results if the diagonals of the adjacency
% matrix are non-zero, for example Aii = degree(i)
% - having an estimation on the upper bound of community numbers "max_rank"
% will significantly increase the performance of the algorithm.
% - P is based on W with zero columns removed.
% - Due to the coordinate descend optimisation scheme, the NMF result is
% initialisation depend. You may have to run the algorithm multiple times if you
% didn't get the desired result from the first time.
% - The main function commDetNMF.m acts as a wrapper to the Tan and Fevotte
% code nmf_kl_mos.m from paper "Automatic Relevance Determination in Nonnegative Matrix
% Factorization", SPARS '09
% - email: ioannis.psorakis@eng.ox.ac.uk for questions/bug reports/comments!
%============

if (nargin < 2)
    max_rank = ceil(size(V,1)/2);
end;

N = size(V,1);

%remove non-connected nodes
active_nodes = sum(V)~=0;
active_node_indices = find(active_nodes);
non_active_nodes = ~active_nodes;
non_active_node_indices = find(non_active_nodes);

if sum(non_active_nodes)>0
   V = V(active_node_indices,active_node_indices);
   fprintf('\nNOTE: the following non-connected nodes are removed from the graphs:\n')
   disp(non_active_node_indices)
end

%now scale the incoming matrix to lie in [0,1] for all elements
X = scale01(V);

%set up the initial W,H matrices
if (nargin < 3)
    W0 = rand(size(X,1),max_rank);
end;
if (nargin < 4)
    H0 = rand(max_rank,size(X,2));
end;

%now run the nmf with shrinkage priors
[W, H] = nmf_kl_mos(X, 100, 'HN', 2, W0, H0, 1, 1, 2);

%allocate using greedy group selection
[a,b] = max(W');

%scale the allocations matrix
P = W ./ repmat(sum(W,2),1,max_rank);
P(:,sum(P)==0)=[];

%now go through the groups finding unique sets
if max_rank>1
m=1;
for n=1:max_rank,
    f=find(b==n);
    if (~isempty(f))
        g{m}=f;
        m=m+1;
    end
end
g = g';

else
    g = (1:size(V,1))';
end

% recover original indices
K = length(g);
number_of_disconnected_nodes = sum(non_active_nodes);
if number_of_disconnected_nodes>0
        
    Pnew = zeros(N,size(P,2));
        
    % fix active
    for k=1:K
       g{k} = active_node_indices(g{k});
    end
    Pnew(active_nodes,:) = P;
    
    % fix inactive
    for k=K+1:K+number_of_disconnected_nodes
       g{k} = non_active_node_indices(k-K);
       Pnew(g{k},:) = 0;
    end    
    
    P = Pnew;
end

end
%--------------
function [W, H, Lk, store] = nmf_kl_mos(V,n_iter,prior,method,varargin)
% Implements the code from Tan and Fevotte "Automatic Relevance Determination in Nonnegative Matrix
% Factorization", SPARS '09

% NMF minimizing KL distance through multiplicative updates
% with specified priors and different methods for model order selection
%
% Input:
%   - V: positive matrix data (F x N)
%   - n_iter: number of iterations
%   - init_W: basis (F x K)
%   - init_H: gains (K x N)
%   - switch_W, switch_H: switches (0 or 1)
%   - b: scale param in beta prior
%
% Alternatively,
% [W, H, cost] = nmf_kl_mos(V, n_iter, prior, method, K)
% uses random initializations.


[F,N] = size(V);

if nargin == 9
    W = varargin{1};
    H = varargin{2};
    switch_W = varargin{3};
    switch_H = varargin{4};
    b0 = varargin{5};
    K = size(W,2);
elseif nargin == 4
    K = varargin{1};
    W = abs(randn(F,K)) + ones(F,K);
    H = abs(randn(K,N)) + ones(K,N);
    switch_W = 1;
    switch_H = 1;
end

a0 = 5;
b0 = 2;
switch_beta = 1;
a = a0*ones(K,1); b = b0*ones(K,1);

% Definitions
store.cost       = zeros(1,n_iter);
store.cost_map   = zeros(1,n_iter);
store.invbeta    = zeros(K,n_iter);

% Initializations
invbeta            = ones(K,1);
V_ap               = W*H;

if strcmp(prior,'HN') && method == 1;
    a_post = N/2 + a -1;
    b_post = sum(H.^2,2)/2 + b;
elseif strcmp(prior,'HN') && method == 2;
    a_post = N/2 + F/2 + a - 1;
    b_post = sum(H.^2,2)/2 + sum(W.^2,1)'/2 + b;
elseif strcmp(prior,'G') && method == 1;
    a_post = N + a - 1;
    b_post = sum(H,2) + b;
elseif strcmp(prior,'G') && method == 2;
    a_post = N + F + a - 1;
    b_post = sum(H,2) + sum(W,1)' + b;
elseif strcmp(prior,'mixed') && method == 2;
    a_post = N/2 + F + a - 1;
    b_post = sum(H.^2,2)/2 + sum(W,1)' + b;
else
    error('no such prior');
end

cost = sum(sum(V.*log(V./V_ap)-V+V_ap));
store.cost(1) = cost;
if method == 1
    store.cost_map(1) = cost + sum(b_post./invbeta) + sum(a_post.*log(invbeta)) +  0.5 * sum(W(:).^2);
elseif method == 2
    store.cost_map(1) = cost + sum(b_post./invbeta) + sum(a_post.*log(invbeta));
else
    error('no such method');
end

store.invbeta(:,1)  = invbeta;

epsilon = sqrt(eps);

%h = waitbar(0,'Please wait...');

%start_time = clock;

for iter = 2:n_iter

    %waitbar(iter/1000,h);

    if switch_H
        if strcmp(prior,'HN') || strcmp(prior,'mixed')
            H = H ./ (W' * ones(F,N) + H./repmat(invbeta ,1,N) + epsilon  ) .* (W' *  (V ./ (W*H)));
        elseif strcmp(prior,'G')
            H = H ./ (W' * ones(F,N) + 1./repmat(invbeta ,1,N) + epsilon  ) .* (W' *  (V ./ (W*H)));
        else
            error('no such method');
        end
    end

    if switch_beta
        if strcmp(prior,'HN') && method == 1;
            b_post = sum(H.^2,2)/2 + b;
        elseif strcmp(prior,'HN') && method == 2;
            b_post = sum(H.^2,2)/2 + sum(W.^2,1)'/2 + b;
        elseif strcmp(prior,'G') && method == 1;
            b_post = sum(H,2) + b;
        elseif strcmp(prior,'G') && method == 2;
            b_post = sum(H,2) + sum(W,1)' + b;
        elseif strcmp(prior,'mixed') && method == 2;
            b_post = sum(H.^2,2)/2 + sum(W,1)' + b;
        else
            error('no such prior');
        end
        invbeta = b_post./a_post;
    end

    if switch_W
        if strcmp(prior,'HN') && method == 1;
            W = W ./ (ones(F,N) * H' + W + epsilon) .* ( (V./(W*H)) * H');
        elseif strcmp(prior,'HN') && method == 2;
            W = W ./ (ones(F,N) * H' + W./repmat(invbeta',F,1) + epsilon) .* ( (V./(W*H)) * H');
        elseif strcmp(prior,'G') && method == 1;
            W = W ./ (ones(F,N) * H' + ones(size(W)) + epsilon) .* ( (V./(W*H)) * H');
        elseif strcmp(prior,'G') && method == 2 || strcmp(prior,'mixed') && method == 2;
            W = W ./ (ones(F,N) * H' + 1./repmat(invbeta',F,1) + epsilon) .* ( (V./(W*H)) * H');
        else
            error('no such prior');
        end
    end

    V_ap = W*H;
    cost = sum(sum(V.*log(V./V_ap)-V+V_ap));

    % Fix scales %%%
    %     if (switch_W && switch_H && 0)
    %         for k = 1:K
    %             scale = sqrt(F / N * sum(H(k,:)) / sum(W(:,k)));
    %             W(:,k) = W(:,k) * scale;
    %             H(k,:) = H(k,:) / scale;
    %         end
    %     end

    store.cost(iter) = cost;

    if method == 1
        store.cost_map(iter) = cost + sum(b_post./invbeta) + sum(a_post.*log(invbeta)) +  0.5 * sum(W(:).^2);
    elseif method == 2
        store.cost_map(iter) = cost + sum(b_post./invbeta) + sum(a_post.*log(invbeta));
    else
        error('no such method');
    end

    store.invbeta(:,iter)  = invbeta;


    if(mod(iter,1e3) == 0)
        disp(['Number of iterations completed = ' num2str(iter)]);
        % disp(['Elapsed time = ' num2str(etime(clock, start_time)) ' seconds.']);
    end

end

Lk = (2*N + 2*(a0-1))/b0;
b = 1./invbeta;
end

%--------------
function y = scale01(x)
  
  y = x - min(min(x)) + eps;
  y = y/max(max(y));
  
return;
end