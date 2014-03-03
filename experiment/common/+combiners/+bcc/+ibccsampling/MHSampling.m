classdef MHSampling < handle
    %MHSAMPLING Independent metropolis-hastings
    %   Sample the independent bayesian classifier combination for a number
    %   of samples
    properties
        nMHSamples = 5000;
        
        combiner;
        
        nClassifiers
        nClasses
        
        lambda
        nu
        
        P
        Alpha
        Pi
        
        T
        C
        nKnownPoints
        
        members = []
    end
    
    methods
        
        function obj = MHSampling(combiner, members)
            obj.combiner = combiner;
            
            obj.nClassifiers = combiner.nClassifiers;
            obj.nClasses = combiner.nClasses;
            obj.lambda = combiner.lambda;
            obj.nu = combiner.nu;
            if nargin > 1
                obj.members = members;
            end
        end
        
        function simulate(obj, nDataPoints, nSamples)
            obj.nMHSamples = nSamples;
            obj.genData(nDataPoints);
            expectedT = obj.sample(obj.C+1, obj.T(1:obj.nKnownPoints));
            figure;
            plot(1:length(expectedT), obj.T - expectedT);
            hold all
            %plot(1:length(expectedT), obj.T);
        end
        
        function ConfSample = sampleConf(obj, Alpha)     
            
            ConfSample = zeros(obj.nClasses, obj.nClasses, obj.nClassifiers);
            
            for k=1:obj.nClassifiers
                %row number corresponds to t=1, t=2
                Conf = gamrnd(Alpha(:, :, k), ones(obj.nClasses, obj.nClasses));
                ConfSample(:, :, k) = Conf ./ (sum(Conf, 2)*ones(1, obj.nClasses));
            end
        end 
        
        function genData(obj, nDataPoints)
            obj.combiner.nClassifiers = obj.nClassifiers;
            
            %nDataPoints = 1;
            obj.nKnownPoints = 0;
            
            obj.nClassifiers = 3;
            obj.nClasses = 2;
            obj.Alpha = zeros(2, 2, 3);
            obj.Alpha(:, :, 1) = [20 3; 3 20];
            obj.Alpha(:, :, 2) = [20 3; 3 20];
            obj.Alpha(:, :, 3) = [20 3; 3 20];
            
            obj.Pi = obj.sampleConf(obj.Alpha);
            
            obj.P = [0.5 0.5];
            
            obj.T = obj.combiner.priorT(obj.P, nDataPoints);
            
            U = rand(obj.nClassifiers, nDataPoints);
            obj.C = zeros(3, nDataPoints);
            obj.C(1, :) = U(1, :) < obj.Pi(obj.T+1, 2, 1)';
            obj.C(2, :) = U(2, :) < obj.Pi(obj.T+1, 2, 2)';
            obj.C(3, :) = U(3, :) < obj.Pi(obj.T+1, 2, 3)';
        end
        
        function expectedT = sample(obj, C, Tknown)
            
            nDataPoints = size(C, 2);
            nKnownPoints = size(Tknown, 2);
            
            acceptedT = zeros(obj.nMHSamples, nDataPoints);
            repetitionT = zeros(obj.nMHSamples, 1);
            
            nAccepted = 0;
            nDistinct = 0;
            nD1 = 0;
            nD2 = 0;
            
            logJointPprev = -inf;
            T = [];
            
            pi1 = zeros(4, obj.nMHSamples);
            a1 = zeros(4, obj.nMHSamples);
            
            PiAcc = [];
            AAcc = [];
            
            logQprev = -inf;
            
            A = [];
            
            
%Code for inspecting the probability values
%             Tnew = obj.T;
%             TnewAsIdx = Tnew + 1;
%             piIdx = sub2ind(size(obj.Pi), ones(obj.nClassifiers, 1)*TnewAsIdx, C, ...
%                 (1:obj.nClassifiers)'*ones(1,nDataPoints));
%             %multiply by P_t_i
%             logpCgivenT = sum(log(obj.Pi(piIdx)), 1);
%             logpCT = sum(logpCgivenT) + sum(log(obj.P(TnewAsIdx)))
%             
%              
%             Tnew = 1 - obj.T;
%             TnewAsIdx = Tnew + 1;
%             piIdx = sub2ind(size(obj.Pi), ones(obj.nClassifiers, 1)*TnewAsIdx, C, ...
%                 (1:obj.nClassifiers)'*ones(1,nDataPoints));
%             %multiply by P_t_i
%             logpCgivenT = sum(log(obj.Pi(piIdx)), 1);
%             logpCT = sum(logpCgivenT) + sum(log(obj.P(TnewAsIdx)))          
%             
            
            AlphaPrev = zeros(obj.nClasses, obj.nClasses, obj.nClassifiers);
            for k=1:obj.nClassifiers
                AlphaPrev(:, :, k) = 1./obj.combiner.lambda();
            end
            alphaSigma = ones(obj.nClasses, obj.nClasses, obj.nClassifiers) * 0.1;
            alphaSigmaExplore = ones(obj.nClasses, obj.nClasses, obj.nClassifiers) * 2;
            exploreFrequency = 4; %once every 4 times
            exploreCount = 1;
            
             PiPrev = zeros(obj.nClasses, obj.nClasses, obj.nClassifiers);
%             for k=1:obj.nClassifiers
%                 PiPrev(:, :, k) = 0.3;
%                 PiPrev(1, 1, k) = 0.7;
%                 PiPrev(2, 2, k) = 0.7;
%             end
%             piSigma = ones(obj.nClasses, obj.nClasses, obj.nClassifiers) * 0.01;
                        
            while nAccepted < obj.nMHSamples
            
                %sample p|nu
                P = obj.combiner.priorP();
                pP = combiners.Ibcc.dirpdf(P(1), obj.nu);
                
                %sample alpha|lambda
                
                %Alpha = obj.Alpha;
                %Alpha = obj.combiner.priorAlpha();
                %should we be fiddling with log alpha instead so that we
                %avoid bouncing off the limit at zero?
                Alpha = -1;
                while ~isempty(find(Alpha <= 0))
                    if exploreCount < exploreFrequency
                        Alpha = normrnd(AlphaPrev, alphaSigma);
                        exploreCount = exploreCount + 1;
                    else
                        Alpha = normrnd(AlphaPrev, alphaSigmaExplore);
                        exploreCount = 1;
                    end
                end
                pAlpha = 1;
                for k=1:obj.nClassifiers
                    pAlpha = pAlpha * prod(prod(exppdf(Alpha(:, :, k), 1./(obj.lambda))));
                end
                
                %sample pi|alpha
                %Pi = obj.Pi;
                Pi = obj.combiner.sampleConf(Alpha);
%                 adj = (alphaSigma-(2*alphaSigma*rand(1,1)));
%                 if adj > 0
%                     adj = adj .* (1-PiPrev);
%                 else
%                     adj = adj .* PiPrev;
%                 end
                
                %this is not quite symmetrical!
%                 Pi = PiPrev + adj;
%                 Pi(1, 2, :) = 1 - Pi(1, 1, :);
%                 Pi(2, 2, :) = 1 - Pi(1, 2, :);
%                 Pi(2, 1, :) = 1 - Pi(2, 2, :);
%                 Pi = -1;
%                 while ~isempty(find(Pi < 0)) || ~isempty(find(Pi > 1))
%                     Pi = normrnd(PiPrev, piSigma);
%                 end
                
                pPi = 1;
                for k=1:obj.nClassifiers
                    pPi = pPi * prod(combiners.Ibcc.dirpdf(Pi(:, :, k), Alpha(:, :, k)), 1);
                end
                
                %sample t|p
                [Tsample, TnewAsIdx, postTi] = obj.combiner.sampleT(C, Pi, P, Tknown+1, obj.members);
                logpT = sum(log(postTi.^(Tsample-1))) + sum(log((1-postTi).^(2-Tsample)));
                
                %TnewAsIdx = [Tknown+1 obj.combiner.priorT(P, nDataPoints-nKnownPoints)+1];
                %logpT = sum(log(P(TnewAsIdx)));
                
                %calculate p(c,t|pi,p)
                Tnew = TnewAsIdx - 1;
                piIdx = sub2ind(size(Pi), ones(obj.nClassifiers, 1)*TnewAsIdx, C, ...
                    (1:obj.nClassifiers)'*ones(1,nDataPoints));
                %multiply by P_t_i
                logpCgivenT = sum(log(Pi(piIdx)), 1);
                logpCT = sum(logpCgivenT) + sum(log(P(TnewAsIdx)));
                
                %the last terms cancel out
                logJointPprop = logpCT   + log(pAlpha);% + log(pPi) + log(pP);
                         
                nAccepted = nAccepted + 1;
                
                %log(q Alpha|AlphaPrev) = log(q AlphaPrev|Alpha) so the
                %terms cancel out
                logQprop = logpT;% + log(qAlpha);% + log(pPi) + log(pP);
                if isnan(logQprop) || isinf(logQprop)
                    display('invalid logQprop');
                end
                %probabilities of producing samples cancel out with terms
                %in the joint probability distribution we are sampling
                %just
                acceptance = logJointPprop - logJointPprev + logQprev - logQprop;
                A = [A acceptance];
                
                if acceptance > 0 || isnan(acceptance)
                    T = Tnew;
                    
                    logJointPprop
                    logQprop
                    
                    logJointPprev = logJointPprop;
                    logQprev = logQprop;
                    AlphaPrev = Alpha;
                    
%                     PiAcc = Pi;
%                     AAcc = Alpha;
                    nDistinct = nDistinct + 1;
                    nD1 = nD1 + 1;
                    acceptedT(nDistinct, :) = T;
                else
                    u = rand(1, 1);%sum(log(rand(nDataPoints, 1)));
                    %if u < exp(acceptance)
                    if log(u) <= acceptance
                        T = Tnew;
                        
                        logJointPprev = logJointPprop;
                        logQprev = logQprop;
                        AlphaPrev = Alpha;
                        
                        
%                         PiAcc = Pi;
%                         AAcc = Alpha;
                        
                        nDistinct = nDistinct + 1;
                        nD2 = nD2 +1;
                        acceptedT(nDistinct, :) = T;                        
                    else
                                                
                    end
                end
                
                if isempty(T)
                    continue;
                end
%                 pi1(1, nAccepted) = PiAcc(1, 1, 1);
%                 pi1(2, nAccepted) = PiAcc(1, 2, 1);
%                 pi1(3, nAccepted) = PiAcc(2, 1, 1);
%                 pi1(4, nAccepted) = PiAcc(2, 2, 1);
%                 
%                 a1(1, nAccepted) = AAcc(1, 1, 1);
%                 a1(2, nAccepted) = AAcc(1, 2, 1);
%                 a1(3, nAccepted) = AAcc(2, 1, 1);
%                 a1(4, nAccepted) = AAcc(2, 2, 1);
                repetitionT(nDistinct) = repetitionT(nDistinct) + 1;                
            end
            
            expectedT = sum(acceptedT .* (repetitionT*ones(1, nDataPoints)), 1) ./ nAccepted;
            
%             figure;
%             plot(1:length(A), A);
%             
%             figure;
%             plot(1:obj.nMHSamples, pi1(1,:), '.');
%             hold all;
%             plot(1:obj.nMHSamples, pi1(2,:), '.');
%             hold all;
%             plot(1:obj.nMHSamples, pi1(3,:), '.');
%             hold all;
%             plot(1:obj.nMHSamples, pi1(4,:), '.');
%             hold off;
            
            display(sprintf('number of distinct samples: %i (%i, %i)', nDistinct, nD1, nD2));
            display(sprintf('out of total samples: %i', nAccepted));
            
        end
        
%         function expectedT = sampleRandomWalk(obj, C, Tknown)
%             
%             nDataPoints = size(C, 2);
%             nKnownPoints = size(Tknown, 2);
%             
%             acceptedT = zeros(obj.nMHSamples, nDataPoints);
%             nAccepted = 0;
%             nDistinct = 0;
%             
%             previousPCT = -inf;
%             T = [];
%             
%             pi1 = zeros(4, obj.nMHSamples);
%             a1 = zeros(4, obj.nMHSamples);
%             
%             PiAcc = [];
%             AAcc = [];
%             
%             P = obj.combiner.priorP();
% 
%             %sample t|p
%             Tsample = obj.priorT(P, nDataPoints-nKnownPoints);
% 
%             Tnew = [Tknown Tsample];
% 
%             %sample alpha|lambda
%             Alpha = obj.combiner.priorAlpha();
% 
%             %sample pi|alpha
%             Pi = obj.combiner.sampleConf(Alpha);            
%             
%             startVector = [P Tnew Alpha Pi];
%             %pdf = @(X) 
%             %proppdf = @(
%             %proprnd
%             smpl = mhsample(startVector, nDataPoints, 'pdf', pdf, ...
%                 'proppdf', proppdf, 'proprnd', proprnd);
%             
%             expectedT = smpl(2, :);
%         end
% 
%         
%         function qCT = q(obj, Pi, Alpha, T, C)
%             
%         end
        
    end    
end

