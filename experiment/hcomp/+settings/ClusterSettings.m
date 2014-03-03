classdef ClusterSettings < handle
%UNTITLED Summary of this class goes here
%   Detailed explanation goes here

    properties       

    %CLUSTERING
        %number of clusters and clustering window length
        seqLength = 20;
        K = 2; 
        clusterData = 'post';
        weightClusterInput = false;
    end
end