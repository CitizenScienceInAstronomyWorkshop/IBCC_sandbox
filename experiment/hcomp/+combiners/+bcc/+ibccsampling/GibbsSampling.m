classdef GibbsSampling < handle 
    properties (Constant)
        % number of iterations - not needed?
        %gibbsIterations = 500;
        %number of accepted samples
        sampleInterval = 1; %1; %3
        burnIn = 40; %20; %40
        gibbsSamples = 300; %100; %200
    end
    properties
        
        combiner;
        
        nClassifiers
        nClasses
        nScores
        
        lambda
        nu
        
        agents
        members = []
        ars
        
        initAlphaFromMeans = true;
        
        convThreshold = 10^-2;
    end
    
    methods
        
        function obj = GibbsSampling(combiner, nClassifiers, members)
            obj.combiner = combiner;
            
            obj.nClassifiers = nClassifiers;
            obj.nClasses = combiner.nClasses;
            obj.nScores = combiner.nScores;
            obj.lambda = combiner.lambda;
            obj.nu = combiner.nu;
            
            obj.agents = combiner.agents;
            if nargin > 2
                obj.members = members;
            end
            
            obj.ars = combiners.bcc.ibccsampling.ARS(obj.nClassifiers, ...
                obj.nClasses, obj.nScores, obj.lambda, obj.nu); 
        end
        
        function [Alpha, P, T, Tknown] = initVariables(obj, Tknown, nSamples)
            %set initial values 
            %start with means of priors or random draws?
            if nargin > 3 && obj.initAlphaFromMeans
                Alpha = zeros(obj.nClasses, obj.nScores, obj.nClassifiers);
                for k=1:obj.nClassifiers
                    Alpha(:, :, k) = 1 ./ obj.lambda;
                end
                P = obj.nu ./ sum(obj.nu); 
            else
                Alpha = obj.combiner.priorAlpha();
                obj.combiner.setAlphaPrior(Alpha, true, obj.combiner.nAgents);
                Alpha = obj.combiner.priorAlpha();
                P = obj.combiner.priorP();
            end
            
            if sum(obj.nu) == 0
                P = [0.5 0.5];
            end            
            
            % no longer necessary as we expect labels in this form now: Tknown = Tknown + 1;
            T = obj.combiner.priorT(P, nSamples - length(find(Tknown~=0))) + 1;
            %T = [Tknown T];            
            %display('more breaking here');
            Tnew = Tknown;
            Tnew(find(Tnew==0)) = T;
            T = Tnew;
        end
        
        function [expectedPost, expectedConf, Alpha] = sample(obj, C, Tknown)
            nSamples = length(Tknown);
                        
            [Alpha0, P, T, Tknown] = obj.initVariables(Tknown, nSamples);
            Alpha = Alpha0;
            currentInterval = obj.sampleInterval;
            nSamplesCollected = 0;
            
            expectedPost = zeros(1, nSamples);
            expectedConf = zeros(obj.nClasses, obj.nScores, obj.nClassifiers);
            %expectedP = zeros(1, obj.nClasses);
            
            i = 0;           
            currentET = 0.5;
%             sumET = 0;

            conv =obj.gibbsSamples;
            while nSamplesCollected < obj.gibbsSamples && i < 5*(obj.gibbsSamples*obj.sampleInterval+obj.burnIn)
                i = i+1;
                display(['gibbs samples collected: ' num2str(nSamplesCollected) ', ' num2str(i)]);

                %Conf conditional is prop. to product of dirichlets
                %p(pi|alpha)p(pi|c,t) = p(pi| alpha + confCounts)
                %confCounts = obj.confusionCounts(T, C, nSamples);
                ET = sparse(obj.nClasses, nSamples);
                for n=1:length(T)
                    %leave the unknown ones as equal probability of each
                    if T(n) == 0
                        %if there is no score for this data point, ignore and set to zero
                        if sum(C{2}==n)==0
                            ET(:, n) = 0;
                        else
                            ET(:, n) = 1 ./ obj.nClasses;
                        end
                    else
                        %with continuous multi-class this goes wrong
                        ET(T(n), n) = 1;
                    end
                end
                confCounts = obj.combiner.voteCounts(C, ET);
                                
                Conf = obj.combiner.sampleConf(Alpha + confCounts, C);
                               
                %uses lambda as prior
                newAlpha = obj.sampleAlpha(Conf, Alpha); 
                if ~isempty(newAlpha)
                    Alpha = newAlpha;
                else
                    %go back to the top of the loop
                    continue;
                end
                
                %need to adapt this so that we include some correct values
                %of T. What about combination on the fly?
                T = obj.combiner.sampleT(C, Conf, P, Tknown, obj.members);
                
                P = obj.sampleP(T); %also uses nu as a prior
                   
                oldET = currentET;                
                
                if i>obj.burnIn && currentInterval >= obj.sampleInterval
                                        
                    if nSamplesCollected==0
                        expectedPost = zeros(1, nSamples);
                    end
                    
                    %T records labels as 1 or 2
                    expectedPost = expectedPost + T - 1;
                    expectedConf = expectedConf + Conf;
                    %expectedP = expectedP + P;
                    currentInterval = 1;
                    nSamplesCollected = nSamplesCollected + 1;
                    %display(sprintf('gibbs sample - %d', nSamplesCollected));
                    
                    currentET = expectedPost ./ nSamplesCollected;
                    if ~isempty(obj.combiner.progressMonitor)
                        obj.combiner.progressMonitor.combinerIteration(...
                            currentET+1, obj.combiner.getId(), i, ...
                            obj.combiner.nKnown); 
                    end
                else
                    expectedPost = expectedPost + T - 1;
                    
                    
                    currentInterval = currentInterval + 1;         
                    if i>obj.burnIn
                        currentET = expectedPost ./ (nSamplesCollected+1);
                    else
                        currentET = expectedPost ./ i;
                    end
                    
                    if ~isempty(obj.combiner.progressMonitor)
                        obj.combiner.progressMonitor.combinerIteration(...
                            currentET+1,...
                            obj.combiner.getId(), i, ...
                            obj.combiner.nKnown); 
                    end
                    if nSamplesCollected > 0
                        expectedPost = expectedPost - T + 1;
                    end
                    
                end               
                
                change = sum(abs(currentET - oldET));
                display(num2str(change));
                if change < obj.convThreshold && conv > i
                    conv = i;
                elseif change >= obj.convThreshold
                    conv = obj.gibbsSamples;
                end
                
                %if no change for 20 samples then assume convergence has
                %been reached
                if conv < i-20
                    break;
                end
                
                if nSamplesCollected==obj.gibbsSamples
                    display('reached the limit');
                end
            end  
            
            nItFile = sprintf('%s%s', obj.combiner.confSave, 'nIts_gibbs');
            if exist(nItFile, 'file');
                prevNIts = dlmread(nItFile);
                prevNIts = [prevNIts i];
            else
                prevNIts = i;
            end
            dlmwrite(nItFile, prevNIts);            
                        
            if nSamplesCollected < obj.gibbsSamples
                display(sprintf('Collected %i.', nSamplesCollected));
            end
            
            %display(sumET./(size(T,2)*20)); %mean difference per data point over the last 20 iterations - was it still  changing significantly?
            expectedPost = expectedPost ./ nSamplesCollected;
            expectedConf = expectedConf ./ nSamplesCollected;
            %expectedP = expectedP ./ nSamplesCollected;
            
            %expectedPost = obj.flipIfNeeded(expectedPost, expectedConf, expectedP, C, Tknown, Alpha0);            
        end  
        
        function EPost = flipIfNeeded(obj, EPost, EPi, EP, C, Tknowns, Alpha)
            [ET, pT] = obj.expectedT(Tknowns, C, EP, EPi);
            logPjoint = sum( log(sum(ET(:, 1:length(Tknowns)) .* pT(:, 1:length(Tknowns)), 1) ), 2)...
                + sum(sum(log(obj.combiner.dirpdf(EPi, Alpha)), 1), 3)...
                + log(obj.combiner.dirpdf(EP, obj.nu));
            
            %flip everything
            EPf = EP;
            EPif = EPi;
            for j=1:obj.nClasses
                
                EPf(j) = EP(mod(j, obj.nClasses)+1);
                EPif(:, j) = EPi(:, mod(j, obj.nClasses)+1);
            end
            
            [ETf, pTf] = obj.expectedT(Tknowns, C, EPf, EPif);
            logPjointf = sum( log(sum(ETf(:, 1:length(Tknowns)) .* pTf(:, 1:length(Tknowns)), 1) ), 2)...
                + sum(sum(log(obj.combiner.dirpdf(EPif, Alpha)), 1), 3)...
                + log(obj.combiner.dirpdf(EPf, obj.nu));
            
            if logPjointf >= logPjoint
                EPost = 1 - EPost;
            end
        end
        
        function [ET, pT] = expectedT(obj, Tknowns, C, lnP, lnPi)
            logJoint = zeros(obj.nClasses, size(C, 2));
            pT = zeros(obj.nClasses, size(C, 2));
            
            [I, J, nonZeroC] = find(C);
            nonZeroIdx = sub2ind(size(C), I, J);
            
            if size(I, 2) > size(I, 1)
                I = I';
            end
            if size(nonZeroC, 2) > size(nonZeroC, 1)
                nonZeroC = nonZeroC';
            end
            for j=1:obj.nClasses
                lnPCT = zeros(1, size(C,2));
                for k=1:obj.nClassifiers
                    lnPCT = lnPCT + lnPi(j, C(k,:), k);
                end

                logJoint(j, :) = lnP(j) + lnPCT;
            end
            pT = exp(logJoint);
            normTerm = ones(obj.nClasses, 1)*sum(pT, 1);
            pT = pT ./ normTerm;
            
            ET = pT;
            ET(:, Tknowns~=0) = 0;
            ET(sub2ind(size(pT), Tknowns(Tknowns~=0), find(Tknowns~=0))) = 1;
        end                
        
        function Alpha = sampleAlpha(obj, Conf, Alpha)
            if obj.combiner.fixedAlpha
                Alpha = obj.combiner.Alpha0;
            else
                Alpha = obj.ars.sampleAlpha(Conf, Alpha);
            end            
        end
        
%         function confCounts = confusionCounts(obj, T, C, nSamples)
%             confusion counts: the observations/predictions produced by 
%             the current model giving the observed confusion counts            
%             confCounts are the counts for each classifier output
%             given each true label                        
%             confCounts = ones(obj.nClasses, obj.nScores, obj.nClassifiers);       
%             
%             T = ones(obj.nClassifiers, 1) * T;
%             
%             %assuming 2 classes
%             for j=1:obj.nClasses
%                 for s=1:obj.nScores
%                     if obj.combiner.scoreSet
%                         confCounts(j, s, :) = sum((C{3}.*(T(C{2})==j)) ==s, 2);                        
%                     else 
%                         confCounts(j, s, :) = sum((C.*(T==j)) ==s, 2);
%                     end
%                 end
%             end
%         end        
        
        function Psample = sampleP(obj, T)
            % p(P|Alpha, conf, T, C) = p(P| T) = p(T|P)p(P)/p(T) prop. to 
            % p(T|P)dirichlet(nu) = P(T=Tsample)dirichlet(nu) =
            % dirichlet(T+nu)
            params = obj.nu;
            for c = 1:obj.nClasses
                params(c) = params(c) + sum(T==c);
            end
                        
            scaleFactors = ones(1, obj.nClasses);
            
            %repeat this a few times so we can see whether it's correct.
            if length(T) > 10
                nPSamples = 10;
            else
                nPSamples = 1;
            end
            
            for i=1:nPSamples
                Psample = gamrnd(params, scaleFactors);
                Psample = Psample ./ sum(Psample);
            end
        end
    end
end