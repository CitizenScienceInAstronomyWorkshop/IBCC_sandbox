classdef logisticRegressionAgent < handle
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
        nTrainingSamples = 0
        
        X %test data
        Ttest %test labels
        
        calibration = 1;   
        bias = 0;
    end
    
    methods
        
        function obj = logisticRegressionAgent(sensors, label, lengthHistory, calibration, bias)
            obj.s = sensors;
            obj.label = label;
            obj.lengthHistory = lengthHistory;
            
            obj.num_active_sensors = 0;
            for s=1:numel(sensors)
                sensor_selected = sensors(1, s);
                if sensor_selected == 1
                    obj.num_active_sensors = obj.num_active_sensors + 1;
                end            
            end
            
            if exist('calibration','var')
                obj.calibration = calibration;
            end
            
            if exist('bias', 'var')
                obj.bias = bias;
            end
        end
        
        function s = getSensors(obj)
            s = obj.s;
            return;
        end   
            
        function train(obj, D, T)
            obj.D = D;
            obj.T = T;
            obj.nTrainingSamples  = size(obj.D,1);
            obj.run_training(D, T);
        end
        
        function run_training(obj, D, T)
            
            T = T(T~=-1);
            D = D(T~=-1, :);            
                        
            obj.W = zeros( obj.num_active_sensors+1, 1 );

            %l = -1; %log likelihood
            delta_w = obj.num_active_sensors; %change in w over each iteration

            N = numel(T);

            %get sub-set of data from the selected sensors
            phi = zeros(N, obj.num_active_sensors+1) ;
            
            %set first column to 1s
            phi(1:N, 1) = ones(N, 1);
            
            sIdx = 1;
            for i=1:length(obj.s)
                sensor_selected = obj.s(1, i);
                if sensor_selected == 1
                    
                    col = D(1:N, i);
                                        
                    phi(1:N, sIdx+1) = col;  % add one as we have inserted a column of zeros
                    sIdx = sIdx + 1;
                end
            end

            error_grad = 1;

            iterations = 0;
            
            while error_grad > 10^-23 && delta_w > 0.001 * numel(obj.s)

                iterations = iterations + 1;
                
                A = phi * obj.W;

                ebx = exp(-A);

                Y = 1.0 ./ (1.0 + ebx);

                error_gradient = transpose(phi) * (Y - T);

                R = zeros(N, N);

                for n=1:N
                    R(n, n) = Y(n, 1) * (1 - Y(n, 1));
                end
                
                H = transpose(phi) * R * phi;
                
                if cond(H) > exp(15)
                    %display(sprintf('break on cond(H) = %i', cond(H) ));
                    break;
                end
                
                delta_w = H\error_gradient;
                
                obj.W = obj.W - delta_w;

                delta_w = abs(delta_w);
                delta_w = sum(delta_w);

                error_gradient = abs(error_gradient);
                error_grad = sum(error_gradient);
                
                %display(sprintf('%i %i %i %i', iterations, error_grad, delta_w, cond(H) ));
                %obj.W

                %l = transpose(T) * log(Y) + (1-transpose(T)) * log(1-Y);
            end
        end
        
        function post = test(obj, X, Ttest)
            obj.init_test(X,Ttest);
            
            N = size(X,2);
            
            post = zeros(1,N);
            
            for n=1:N
                post(n) = classify(n);
            end
        end
        
        function init_test(obj, X, Ttest)
            obj.X = X;
            obj.Ttest = Ttest;
        end
        
        function [posterior, a] = classify(obj, n)
            
            d = obj.X(n, :);
            
            phi = ones(1, obj.num_active_sensors + 1);
            
            sIdx = 1;
            for i=1:length(obj.s)
                %activation = obj.s{1, i};
                if obj.s(1, i) == 1
                    col = d(1, i);
                    phi(1, sIdx+1) = col;  % add one as we have inserted a column of zeros
                    sIdx = sIdx + 1;
                end
            end
            a = phi * obj.W;
            ebx = exp(-a);
            posterior = 1.0 / (1.0 + ebx);
                        
            if obj.bias>0
                display('applying bias to agents before calibrating');
            end
            posterior = posterior + obj.bias;
            if posterior > 1
                posterior = 1;
            elseif posterior < 0
                posterior = 0;
            end 
            
            strength = posterior -0.5;
            if strength>0
                l = 1;
            else
                l = -1;
            end
            posterior = 0.5 + l* abs(strength).^obj.calibration;
            if posterior > 1
                posterior = 1;
            elseif posterior < 0
                posterior = 0;
            end                
        end
        
        function retrain(obj, n)
            %number of points we need from the original training set
            deficit = -n+obj.lengthHistory;
                        
            if deficit > 0
                
                if deficit > obj.nTrainingSamples
                    start = 1;
                else 
                    start = obj.nTrainingSamples - deficit;
                end
                
                training_portion = obj.D(start:obj.nTrainingSamples, :);
                test_portion = obj.X(1:n, :);
                
                D = [training_portion; test_portion];
                
                training_labels = obj.T(start:obj.nTrainingSamples, :);
                test_labels = obj.Ttest(1:n, :);
                
                T = [training_labels; test_labels];
                          
            else
                D = obj.X(n-obj.lengthHistory+1:n, :);
                T = obj.Ttest(n-obj.lengthHistory+1:n, :);
            end
            

            obj.run_training(D, T);

        end
    end
    
end

