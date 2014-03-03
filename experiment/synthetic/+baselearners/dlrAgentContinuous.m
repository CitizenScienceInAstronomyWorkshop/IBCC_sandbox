classdef dlrAgent_continuous < handle
    %DLR agent that does not retrain the whole data set but continues to
    %update the existing weights.
    
    properties
        W %weights

        s %binary vector listing which sensors to use
        label %index or name for identifying the agent
        D %training data saved for retraining
        T %labels for saved training data
        num_active_sensors %number of active sensors        
        dlr
        
        X %test data
        Ttest %test labels        
    end
    
    methods
        
        function obj = dlrAgent(sensors, label)
            obj.s = sensors;
            obj.label = label;

            obj.num_active_sensors = sum(sensors);
            
            obj.dlr = dlrpadapter();
        end
            
        function train(obj, D, T)
            obj.D = D;
            obj.T = T;
            
            obj.dlr.dlrpadapt(obj.D(obj.s, :), obj.T);
        end
        
        function [posterior, a] = classify(obj, d)
            d = obj.X(n, :);
            d = d(obj.s, :);
            
            %process with unknown label
            os = obj.dlr.process_data_point(d, -1);
            posterior = os.y(length(os.y));
            a = os.a(length(os.a));
        end
        
        function init_test(obj, X, Ttest)
            obj.X = X;
            obj.Ttest = Ttest;
                      
        end
        
        function retrain(obj, n)
            d = obj.X(n, :);
            
            d = d(obj.s, :);
            
            %process with known label
            obj.dlr.process_data_point(d, t);           
        end
    end
    
end

