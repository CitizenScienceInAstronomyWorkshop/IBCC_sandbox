classdef Weights < handle
    %WEIGHTS Summary of this class goes here
    %   Detailed explanation goes here
        
    methods (Abstract)
        W = calculateWeights(obj, scores, post, T, combinedPost, n);
    end
    
end

