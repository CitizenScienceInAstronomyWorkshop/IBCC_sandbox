classdef ARS < handle
    %adaptive rejection sampling for sampling the alpha conditional in
    %Ibcc.
    %This still causes problems, i.e. doesn't work very well.
    properties
        nStupidAlpha = 0
        nSensibleAlpha = 0
        
        nClassifiers
        nClasses
        nScores
        
        lambda
        nu
        
        maxAlpha = 100; %artificially limit alpha to a sensible amount
        minAlpha = 1e-07
        
        Z
        X
        hX
        gradHX
        uZ
        zCD
    end
    
    methods
        
        function obj = ARS(nClassifiers, nClasses, nScores, lambda, nu)
            obj.nClassifiers = nClassifiers;
            obj.nClasses = nClasses;
            obj.nScores = nScores;
            obj.lambda = lambda;
            obj.nu = nu;
        end
                
        function gradHX = gradHx(obj, X, sumAlphaOthers, logPi, l)            
            %AlphaOthers are the other parameters in this row of pi 
                        
            if X < obj.minAlpha
                digamX = inf;
            else
                digamX = psi(X);
            end
                        
            D = psi(sumAlphaOthers + X);
            gradHX = D - l + logPi - digamX;

            if isinf(gradHX)
                display('gradient of H(X) is infinite');
            end
        end
                
        function [gAlpha, hAlpha] = g(obj, alpha, sumAlphaOthers, pi, l)
            hAlpha = obj.h(alpha, sumAlphaOthers, pi, l);
            gAlpha = exp(hAlpha);
        end
        
        function hAlpha = h(obj, alpha, sumAlphaOthers, pi, l)
            
            gammaSumTerm = log(gamma(sumAlphaOthers + alpha));            
            
            hAlpha = gammaSumTerm - l.*alpha + log(pi.^(alpha-1)) - log(gamma(alpha));          
            
            if isinf(hAlpha)
                display('infinite h(alpha)');
            end    
            if isnan(hAlpha)
                display('hAlpha is not a number');
                display(pi);
                display(alpha);
                display(sumAlphaOthers);
            end
        end                  
        
        function [Alpha] = sampleAlphaElement(obj, k, j, c, l, pi, Alpha0)
            
            if isnan(pi)
                Alpha = 0;
                obj.hX = 0;
                obj.gradHX = 0;
                return;
            end
            
            nGrid = 5;
            lc = l(c);
            %logLc = log(lc);
            
            logPic = log(pi(c));
                                      
            maxX = 20;
            if logPic == -inf
                nGrid = 1;
                obj.X = 1;
                display(sprintf('log pi(%d) for classifier %d and true label %d is -inf', c, k, j));                
            elseif logPic == 0
                display(sprintf('log pi(%d) for classifier %d and true label %d is 0', c, k, j));
                obj.X = (1:nGrid) .* maxX ./ nGrid;
            else
                obj.X = (1:nGrid) .* maxX ./ nGrid;
            end
            
            lGrid = lc * ones(1, nGrid);
            picGrid = pi(c) * ones(1, nGrid);
            logPicGrid = logPic * ones(1, nGrid);
              
            AlphaOthers = Alpha0(j, :, k);
            if c == 1
                AlphaOthers = AlphaOthers(c+1:obj.nScores);
            elseif c==obj.nClasses
                AlphaOthers = AlphaOthers(1:c-1);
            else
                AlphaOthers = [AlphaOthers(1:c-1) AlphaOthers(c+1:obj.nScores)];
            end    
            
            sumAlphaOthers = sum(AlphaOthers);
            
            %The first point X at least
            %should have a positive gradient in the log space
            obj.gradHX = obj.gradHx(obj.X, sumAlphaOthers, logPicGrid, lGrid);

            %make sure at least the last grid point is on the downward slope
            i = 1;
            lastGrad = obj.gradHX(nGrid);
            adj = lastGrad + 1;
            while lastGrad > 0 && obj.X(length(obj.X)) < obj.maxAlpha             
                obj.X = obj.X * adj;
                lastGrad = obj.gradHx(obj.X(nGrid), sumAlphaOthers, logPic, lc);
                i = i+1;
            end
            if i > 1
                obj.gradHX = obj.gradHx(obj.X, sumAlphaOthers, logPicGrid, lGrid);
            end
            
            %while the first grid point still has some components with negative gradient
            maxIt = 100;
            i = 1;               
            while obj.gradHX(1) < 0 && i<maxIt && obj.X(1) > 0.0001 %&& gradHX(1) > -inf
                obj.X(1) = obj.X(1) / 2; 
                obj.gradHX(1) = obj.gradHx(obj.X(1), sumAlphaOthers, logPic, lc);
                i = i+1;
                
                while isinf(obj.gradHX(1))
                    obj.X(1) = obj.X(1) * 2;
                    obj.gradHX(1) = obj.gradHx(obj.X(1), sumAlphaOthers, logPic, lc);
                end
            end
                        
            obj.hX = obj.h(obj.X, sumAlphaOthers, picGrid, lGrid);            
            
            %drop Xs that produce neg inf apart from first one
%             negInfIdx = find(isinf(hX));
%             if ~isempty(negInfIdx)
%                 
%                 newEndIdx = negInfIdx(1) - 1;
%                 
%                 X = X(1:newEndIdx);
%                 hX = hX(1:newEndIdx);
%                 gradHX = gradHX(1:newEndIdx);
%                 nGrid = newEndIdx;
%                 
%                 display(['reduced the grid points due to inf hX to ' num2str(X)]);
%             end
            
            %drop Xs that produce Nan
            nanIdx = find(isnan(obj.hX));
            if ~isempty(nanIdx)
                
                newEndIdx = nanIdx(1) - 1;                
                
                obj.X = obj.X(1:newEndIdx);
                obj.hX = obj.hX(1:newEndIdx);
                obj.gradHX = obj.gradHX(1:newEndIdx);
                nGrid = newEndIdx;
                
                display(['reduced the grid points due to hX=NaN to ' num2str(obj.X)]);
            end
            
            accepted = false;
            
            i = 1;
            
            while ~accepted && i < 100
                [Alpha, idx] = obj.sampleEnvelope();
                [accepted, hAlpha] = obj.acceptSample(Alpha, sumAlphaOthers, pi(c), lc);
                i = i+1;
                if ~accepted && Alpha > 0 && ~isnan(Alpha) && Alpha < obj.maxAlpha
                    if sum(Alpha > obj.X(idx)) > 0
                        xIdx = idx;
                    else
                        xIdx = idx-1;
                    end
                    
                    gradAlpha = obj.gradHx(Alpha, sumAlphaOthers, logPic, lc);
                    
                    if isinf(Alpha)
                        display('got an infinite alpha sample');
                    elseif isinf(exp(Alpha * gradAlpha))
                        display('got an alpha value that was too large to be useable');
                    elseif isinf(gradAlpha)
                        display('got a value of alpha gradient that is inf');
                    else
                        if xIdx == nGrid
                            obj.X = [obj.X Alpha];
                            obj.hX = [obj.X hAlpha];
                            obj.gradHX = [obj.gradHX gradAlpha];
                        elseif xIdx == 0
                            obj.X = [Alpha obj.X];
                            obj.hX = [hAlpha obj.hX];
                            obj.gradHX = [gradAlpha obj.gradHX];                        
                        else
                            obj.X = [obj.X(1:xIdx) Alpha obj.X(xIdx+1:nGrid)];
                            obj.hX = [obj.hX(1:xIdx) hAlpha obj.hX(xIdx+1:nGrid)];
                            obj.gradHX = [obj.gradHX(1:xIdx) gradAlpha obj.gradHX(xIdx+1:nGrid)];
                        end
                        
                        nGrid = nGrid + 1;
                    end                    
                end
            end
            
            if i >= 100 && ~accepted
                Alpha = -1;
                display('rejecting this gibbs iteration totally');
            end
        end
                
        function AlphaSample = sampleAlpha(obj, Conf, Alpha0)                 
            
            %can't easily apply rejection sampling with a multivariate
            %alpha, so we continue the gibbs procedure and sample each
            %element of Alpha in turn using ARS.
                                    
            %start with the initial value for alpha (sample from previous
            %gibbs iteration or from the priors) and update it with new samples. 
            AlphaSample = Alpha0;%zeros(obj.nClasses, obj.nClasses, obj.nClassifiers);
                       
            for k=1:obj.nClassifiers
                for j=1:obj.nClasses
                % We perform rejection sampling using function g(A).
                %First, select a set of points X spread around the peak of
                %lnf(A).                
                    l = obj.lambda(j, :);
                    pi = Conf(j, :, k);                
                    for c=1:obj.nScores
                        alpha = obj.sampleAlphaElement(k, j, c, l, pi, AlphaSample);
                        if alpha <= 0
                            AlphaSample = [];
                            return;
                        end
                        AlphaSample(j, c, k) = alpha;
                    end
                end                
            end         
        end     

        function intersects(obj)
            %calculate Z, the intersections of each segment of the envelope
            nX = size(obj.X, 2);
            
            hXPlusOne = obj.hX(2:nX);
            hXtemp = obj.hX(:, 1:nX-1);
            
            XPlusOne = obj.X(:, 2:nX);
            Xtemp = obj.X(:, 1:nX-1);

            gradHXPlusOne = obj.gradHX(:, 2:nX);
            gradHXtemp = obj.gradHX(:, 1:nX-1);
            
            constDiff = hXPlusOne - XPlusOne.*gradHXPlusOne - hXtemp + Xtemp.*gradHXtemp;
            gradDiff = gradHXtemp - gradHXPlusOne;
            obj.Z = constDiff ./ gradDiff;
            if sum(isinf(obj.Z)) > 0
                obj.Z(isinf(hXPlusOne)) = obj.X(isinf(hXPlusOne));
            end
            
            %add the last value of Z back on - we don't calculate this as
            %alpha is unbounded so the highest value of Z is infinity.
            obj.Z = [obj.Z inf];
            
            if sum(isinf(obj.Z)) > 1
                display('inf z');
            end
            
        end     
        
        function uAlpha = lnEnvelope(obj, Alpha, ZIdx)
            %find out which section applies to Alpha?
            Xj = obj.X(ZIdx);
            hXj = obj.hX(ZIdx);
            dist = Alpha - Xj;
            uAlpha = hXj + dist.*obj.gradHX(ZIdx);
            
            if isnan(obj.uZ)
                display('Alpha candidate sample is not a number');
                obj.gradHX(ZIdx)
            end
        end        
        
        function zCumDensity(obj)
            %cumulative density for the segments of the envelope between z-1 and z
            %Z0 = 0 because alpha can't be lower than 0
            
            %add on the left hand point, z0 and forget the right hand
            %point, which is equal to infinity
            
            nZ = size(obj.Z, 2);
            
            if nZ < 2
                obj.zCD = 1;
                return;
            end
            
            Ztemp = [0 obj.Z(:, 1:nZ-1)];
            gradHXtemp = obj.gradHX(:, 1:size(obj.gradHX, 2)-1);
            uZright = obj.uZ(:, 1:nZ-1);
            
            Zleft = Ztemp(:, 1:nZ-1);
            Zright = Ztemp(:, 2:nZ);
            
            %constant
            b = uZright - gradHXtemp.*Zright;
            
            %Diff = exp(gradHX.*Zright) - exp(gradHX.*Zleft);
            Diff = exp(gradHXtemp.*Zleft) - exp(gradHXtemp.*Zright);
            zD = exp(b)./-gradHXtemp .* Diff;
           
            %need to add the remaining density for the long tail
            obj.zCD = cumsum(zD);
            obj.zCD = obj.zCD ./ obj.zCD(nZ-1);
            obj.zCD = [obj.zCD 1];
        end        
        
        function idx = selectZ(obj)
            %Return index of a random Z value; weight by
            %density when selecting
            U = rand(1, 1) * ones(1, size(obj.zCD, 2));
            
            possIdxs = find((obj.zCD - U)>0);
            idx = possIdxs(1);
        end
        
        function [Alpha, idx, Z] = sampleEnvelope(obj)
            %sample from the envelope distribution
            
            if isinf(obj.hX)
                display('we have an infinite hX, possibly caused by pi=0');
                %set alpha to minus 1 - reject it!
                Alpha = -1;
                idx = 1;
                Z = [];
                return
            end
            
            %find intersects
            obj.intersects();
            obj.uZ = obj.lnEnvelope(obj.Z, 1:size(obj.Z, 2));          
            obj.zCumDensity();
            idx = obj.selectZ();

            %sample from the chosen Z    
            U = rand(1, 1);
            
            grad = obj.gradHX(idx);
            if idx == 1
                eZMinus1 = 1;
            else
                eZMinus1 = exp(grad*obj.Z(idx-1));
            end
            eZ = exp(grad*obj.Z(idx));
            
            if isinf(eZ) && eZ > 0 % grad > 0 and Zidx is inf. eZ = +inf
                display('pos gradient on the last segment. Not sure this will work');
                %replace with e(maxAlpha) as this is the largest we can
                %take anyway. 
                eZ = exp(grad*obj.maxAlpha);
            end
            
            if isinf(grad)
                if idx ==1
                    Alpha = U * obj.X(idx)-obj.Z(idx-1);
                else
                    Alpha = U * obj.X(idx);
                end
                obj.nStupidAlpha = obj.nStupidAlpha + 1;
            elseif eZMinus1 == 0 % && isinf(eZ); don't check this as it's always true if the first condition is true
                sub = -100;
                
                factor = grad*obj.Z(idx-1) - sub;
                
                eZMinus1 = exp(sub);
                Alpha = (factor + log(U*(eZ - eZMinus1) + eZMinus1))/grad;
                
                obj.nStupidAlpha = obj.nStupidAlpha + 1;
            else
                obj.nSensibleAlpha = obj.nSensibleAlpha + 1;
                
                Alpha = log(U*(eZ - eZMinus1) + eZMinus1)/grad;
            end
            if isinf(Alpha)
                display('yarg. inf alpha just produced.');
            elseif isnan(Alpha)
                display('not a number alpha sample!');
            elseif Alpha < 0
                display('somehow have a neg alpha');
            end
        end
        
        function zIdx = findSegment(obj, Alpha)
            nSamples = size(Alpha, 2);
            nZ = size(obj.Z, 2);
            
            zIdx = zeros(1, nSamples);
            
            for i=1:nSamples
            
                A = Alpha(i) * ones(1, nZ);

                diff = obj.Z - A;

                %Z has index 1 to nZ
                % we also have Z0 = -inf and Z_nZ+1 = +inf

                Zgreater = find(diff>0); %find the idxs of Z > X    
                %only need the adjacent value of Z
                if isempty(Zgreater)
                    display('arghh!!!! A is greater than/equal to infinity?');
                    display(A);
                    zIdx(1, i) = size(obj.Z, 2);
                else
                    zIdx(1, i) = Zgreater(1);
                end
            end
        end        
        
        function [sk, uk] = envelope(obj, Alpha)
            zIdx = obj.findSegment(Alpha);
            uk = obj.lnEnvelope(Alpha, zIdx);
            sk = exp(uk);
        end
        
        function [accept, hAlpha] = acceptSample(obj, Alpha, AlphaOthers, pi, l)
            % calculate g(Alpha) and env(Alpha) and decide whether to accept
                                  
            [gAlpha, hAlpha] = obj.g(Alpha, AlphaOthers, pi, l);
            
            if Alpha < obj.minAlpha || Alpha > obj.maxAlpha
                accept = false;
                return;
            end
            
            if isinf(hAlpha)
                %just accept as there is no alternative
                accept = true;
                return
            end
            
            if isnan(hAlpha)
                accept = false;
                return
            end
            
            [eAlpha, uAlpha] = obj.envelope(Alpha);
            
            U = rand(1, 1);
            
            %e^h / e^u = e (h-u); lnU = h-u
            if hAlpha - uAlpha > log(U)
                accept = true;
            else
                accept = false;
            end
        end
    end
end