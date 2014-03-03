classdef dlrpadapter < handle
    
    properties
        os
        t
    end
    
    methods
        
        function obj = dlrpadapter()
            t = 0;
        end
        
        function os = process_data_point(obj, x, label)
            
            obj.t = obj.t+1;
            t = obj.t;
            
            os = obj.os;    
            
            if(length(os.z) < t)
                %resize
                os.y = [os.y; 0];
                os.yolds = [os.yolds; 0];
                os.a = [os.a; 0];
                os.z = [os.z; 0];
                os.pc = [os.pc; 0];
                os.ar = [os.ar; 0];
                os.w = [os.w; zeros(1, D+1)];
                os.pnoises = [os.pnoises; 0];                
            end
            
            if (rem(t,100) == 0)
                fprintf('.'); %print re-assuring dots
            end;
            if (t > BurnIn)
              temp = mean(abs(os.yolds(1:t) - os.z(1:t)) > 1 - epsilon);
              denom = 0;
              for i=1:numbins
                denom = denom + ...
                  (1/2+(1/2)*(1/numbins)) * ...
                  mean((min(os.yolds(1:t),1-os.yolds(1:t)) < i*pnoise/(2*numbins)) .* ...
                       (min(os.yolds(1:t),1-os.yolds(1:t)) >= (i-1)*pnoise/(2*numbins)));
              end
              if (denom > 0)
                  pnoise = temp / denom;
              end;
            end
            pnoise = max(pnoise,0.01); %have a low lower bound
            os.pnoises(t) = pnoise;
            
            Sp = S + q*I; % current adjusted sigma
            [yold,a,s] = lat2out(w,Sp,x); % get output from latent space
            y = (1-pnoise)*yold + pnoise/2;

            os.yolds(t) = yold;
            os.a(t) = a;

            if (label == -1) % no label is known
              z = y;
              dq = y*(1-y); 
            else
              z = label;
              dq = 0;
            end;
            u = y*(1-y);
            Snew = Sp -( (y*(1-y))/(1 + y*(1-y)*s) )*(Sp*x)*(Sp*x)'; 
            %wnew = w + Snew*x*(z-y); % Newton step
            wnew = w + (Sp/(1+y*(1-y)*s))*x*(z-y); % EKF step
            S = Snew;
            w = wnew;
            [yup,aup,sup] = lat2out(w,Snew,x);
            yup = (1-pnoise)*yup + pnoise/2;
            uup = yup*(1-yup);
            v = v + (uup-u);
            if  (Adapt == 1)
              q = max(uup - u, 0) + dq;
            end;

            os.z(t) = z;
            os.y(t) = y;
            os.pc(t) = y>0.5;
            os.ar(t) = mean(trace(S));
            os.w(t,:) = w;  
            
            obj.os = os;
        
    
        function os = train(obj, X, T)
            %function os = dlrpadapt(X, T);
            % X inputs, T targets : 0 or 1 {classes] or -1 [no label known]

            % os : output structure with fields:
            %   y       : pr(class)
            %   yolds   : pr(class) assuming no noise
            %   a       : activation assuming no noise
            %   z       : labels used for training
            %   w       : weights
            %   pc      : predicted class labels
            %   ar      : effective adaption rate of the process
            %   pnoises : estimated label noise

            % (c) Stephen Roberts, Duncan Lowne, Roman Garnett (2007-9)

            Adapt = 1; %sets for adaptation of the parametric weights or not [1/0]
            BurnIn = 100; %used for the evaluation of the estimated label noise level

            [L,D] = size(X);

            if  (L<D)
              X = X';
              [L,D] = size(X);
            end

            pnoise = 0.1; %this is the prior

            if (length(T) ~= L)
              error('X,T must be same length');
            end

            % now read in everything and can start processing

            % normalise the entire data set and augment with column of ones (for bias weights)
            f = find(X == 0);
            X = normalis(X,X);
            X(f) = 0;
            X = [X , ones(L,1)];
            dim = D+1; % dimensionality of weights

            % first need to set up the priors
            I = eye(dim);
            Ir = diag([ones(1,dim-1),0]); % no noise on bias term latent updates if
                                          % we use this
            S = I;
            S2 = S*100;
            w = zeros(dim,1); % last w is the bias
            w2 = w;
            q = 0;
            q2 = q;
            v=0;
            v2 = v;
            dq = 0;
            dq2 = dq;

            %now pre-allocate for speed
            os.y = zeros(L,1);
            os.yolds = zeros(L,1);
            os.a = zeros(L, 1);
            os.z = zeros(L,1);
            os.pc = zeros(L,1);
            os.ar = zeros(L,1);
            os.w = zeros(L, D+1);
            os.pnoises = zeros(L,1);

            obj.os = os;

            numbins = 6;
            epsilon = 0.1;

            for t=1:L % scroll through the data
                x = X(t,:)'; % current input
                obj.process_data_point(x, T(t));
            end;
            
            os = obj.os;
            
            fprintf('\n');
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        function [y,a,s] = lat2out(w,S,x)
            a = w'*x; % latent activation
            s = x'*S*x;
            y = moderate(a,s);
        end
    end
end
   
