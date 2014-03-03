Implements the community detection methodology presented in 
"Overlapping Community Detection using Nonnegative Matrix Factorization", Physical Review E 83, 066114 (2011)
by Ioannis Psorakis, Stephen Roberts, Mark Ebden and Ben Sheldon.
University of Oxford, United Kingdom.

Inputs :
V: NxN adjacency matrix
max_rank: a prior estimation on the maximum number of communities
(optional).

Outputs :

P: NxK matrix where each element Pik represents the (normalised)
soft membership score of node-i to each community-k. P captures the OVERLAPPING
community structure.

g: groups cell array. Each element g{k} contains the node indices of
community k. Each node i is placed in the community g{k} for which it has
the maximum Pnk. Therefore g describes the NON-OVERLAPPING community
structure.

W,H are the NMF factors (normalised from 0 to 1)

NOTES:
------
- In many cases the algorithm gives better results if the diagonals of the adjacency
matrix are non-zero, for example Aii = degree(i)
- having an estimation on the upper bound of community numbers "max_rank"
will significantly increase the performance of the algorithm.
- P is based on W with zero columns removed.
- Due to the coordinate descend optimisation scheme, the NMF result is
initialisation depend. You may have to run the algorithm multiple times if you
didn't get the desired result from the first time.
- The main function commDetNMF.m acts as a wrapper to the Tan and Fevotte
code nmf_kl_mos.m from paper "Automatic Relevance Determination in Nonnegative Matrix
Factorization", SPARS '09
- email: ioannis.psorakis@eng.ox.ac.uk for questions/bug reports/comments!
============