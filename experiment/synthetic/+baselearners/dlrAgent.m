classdef dlrAgent < handle
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        W %weights

        s %binary vector listing which sensors to use
        label %index or name for identifying the agent
        D %training data saved for retraining
        T %labels for saved training data
        num_active_sensors %number of active sensors
        lengthHistory %how many past data points do we use when retraining
        
        os % results of the dlr
        
        nTrainingSamples %size of initial training set, i.e. data before
        %we start looking at predictions
        
        padapt %use dlrpadapt?
    end
    
    methods
        
        function obj = dlrAgent(sensors, label, padapt)
            obj.s = sensors;
            obj.label = label;
            obj.num_active_sensors = sum(sensors);
            obj.padapt = padapt;
        end
            
        function train(obj, D, T)
            %store it for now as we run training and testing in one step
            obj.D = D; %input data
            obj.T = T; %target labels
            obj.nTrainingSamples = length(obj.T);
        end
        
        function init_test(obj, X, Ttest)
            obj.D = [obj.D; X]; %adds test data to training data
            obj.T = [obj.T; Ttest]; %adds test label to training label
            
            obj.D = obj.D(:, find(obj.s));
            
            if obj.padapt
                obj.os = dlrpadapt(obj.D, obj.T);
            else
                obj.os = dlr(obj.D, obj.T);
            end
            
        end
        
        function [posterior, a] = classify(obj, n)
            
            posterior = obj.os.y(n + obj.nTrainingSamples);
            a = obj.os.a(n + obj.nTrainingSamples);
            obj.W = transpose(obj.os.w(n, :));
        end
        
        function retrain(obj, n)
            %do nothing as we already trained at the start
        end
    end
    
end

