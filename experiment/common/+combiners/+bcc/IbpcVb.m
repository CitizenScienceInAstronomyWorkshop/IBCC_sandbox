classdef  IbpcVb < combiners.bcc.Ibcc
    %Independent Bayesian Probability Combination: combines
    %probabilisitic outputs from classifiers (rather than discrete
    %classifications

    properties
        maxIt = 500;
        convIt = 10;
        
        sampleAlpha = false;
        ars
        
        %Extra Hyperhyperparameters
        Mu_lambda; %prior mean for the output of a classifier given a class
        Tau_lambda; %precision in the prior over that mean
        
        Theta_gamma; %scale of precision over the mean a: inverse of how much the mean a-value can vary
        Theta_alpha; %shape of belief in precision of classifier outputs
        Theta_beta; %scale of precision of classifier outputs
    end
 
    methods (Static)
        function sl = shortLabel
            sl = 'Ibpc';
        end   
    end
    
    methods       
        function obj = IbpcVb(nAgents, K, targets, agents, nClasses, nScores, ...
                lambdaMag, lambdaSym, nu)
            obj@combiners.bcc.Ibcc(nAgents, K, targets, agents, nClasses, nScores);

            obj.label = 'I.B.P.C.'; 
            
            obj.ars = combiners.bcc.ibccsampling.ARS(nAgents, ...
                obj.nClasses, obj.nScores, obj.lambda, obj.nu); 
            
            obj.nonDet = false;
            
            if nClasses>2
                display('Not suitable for multiclass yet');
            end
            
            obj.Mu_lambda = ones(nClasses, nAgents);
            %a1 = obj.prepareA(lambdaMag.*(1-lambdaSym));
            %a0 = obj.prepareA(lambdaMag.*lambdaSym./(nClasses-1));
            a0 = -20;
            a1 = 20;
            obj.Mu_lambda(1,:) = a0;
            obj.Mu_lambda(2,:) = a1;
            
            obj.Tau_lambda = ones(nClasses, nAgents);
            obj.Tau_lambda(1,:) = 1 ./ (abs(obj.Mu_lambda(1, :)).^0.5);
            obj.Tau_lambda(2,:) = 1 ./ (abs(obj.Mu_lambda(2, :)).^0.5);
            
            obj.Theta_gamma = ones(nClasses, nAgents);
            obj.Theta_gamma(1,:) = 5;
            obj.Theta_gamma(2,:) = 5;
            
            obj.Theta_alpha = ones(nClasses, nAgents);
            obj.Theta_alpha(1,:) = 5;
            obj.Theta_alpha(2,:) = 5;
            
            obj.Theta_beta = ones(nClasses, nAgents);
            obj.Theta_beta(1,:) = 25;
            obj.Theta_beta(2,:) = 25;        
        end
        
        function A = prepareA(obj, post)
            %logodds values from posteriors
            A = log(post) - log(1-post);
        end
        
        function [Lambda, Gamma, Alpha, Beta, Mu_a, Tau_a, P, T, Tknown] = initVariables(obj, Tknown, A)
                     
            Lambda = obj.Mu_lambda;
            Gamma = obj.Theta_gamma;
            Alpha = obj.Theta_alpha;
            Beta = obj.Theta_beta;
            
            P = obj.nu ./ sum(obj.nu); 
            if sum(obj.nu) == 0
                P = [0.5 0.5];
            end            
            
            Tknown = Tknown + 1;
            T = zeros( obj.nClasses, length(Tknown));
            for n=1:length(Tknown)
                %leave the unknown ones as equal probability of each
                if Tknown(n) == 0
                    T( :, n) = 1 ./ obj.nClasses;
                else
                    T( Tknown(n), n) = 1;
                end
            end
            
            [Mu_a, Tau_a] = obj.expectedPi(A(:, Tknown>0), T(:, Tknown>0), Lambda, Gamma, Alpha, Beta);
        end        
                
        function combinedPost = combineCluster(obj, post, T, members)          
            obj.nKnown = round(obj.nKnown);
            Tknown = T;
            
            A = obj.prepareA(post);
                        
            converged = false;
            pT = 0;
            
            nIt = 0;
            cIt = 0;
            
            %change T into the matrix format for multi-class problems
            [Lambda, Gamma, Alpha, Beta, Mu_a, Tau_a, P, T, Tknowns] = ...
                obj.initVariables(Tknown, A);
            while ~converged
                oldPT = pT;
                
                [T, pT] = obj.expectedT(Tknowns, A, P, Mu_a, Tau_a);
                
                pT = prod(sum(T .* pT, 1), 2);
                                
                nIt = nIt + 1;
                if pT == oldPT
                    cIt = cIt + 1;
                else
                    if cIt > 0
                        display(sprintf('reset cIt: old p(T)=%f, new p(T)=%f', oldPT, pT));
                    end
                    cIt = 0;
                end
                if cIt > obj.convIt  || nIt > obj.maxIt
                    converged = true;
                    display(sprintf('IbPcVB iterations: %d, iterations while converged %d, ln(pt)=%f', nIt, cIt, log(pT)));
                end
                
                %calc expected values
                Lambda = obj.expectedLambda(Mu_a, Tau_a, Gamma);
                Gamma = obj.expectedGamma(Mu_a, Tau_a, Lambda);
                %this will perhaps need to be repeated in an additional
                %inner loop until beta and alpha converge.
                Beta = obj.expectedBeta(Alpha, Tau_a);
                Alpha = obj.expectedAlpha(Beta, Tau_a);
                [Mu_a, Tau_a] = obj.expectedPi(A, T, Lambda, Gamma, Alpha, Beta); %this is slow
                P = obj.expectedP(T);
                                
                if pT < oldPT 
                    display(sprintf('IbccVb: %d, %f', nIt, log(pT)));
                end
            end
            
            %need to change this for >2 classes as it only gives
            %probability of T==2
            combinedPost = T(2, :);
            
            %save the expected confusion distributions
            dlmwrite(sprintf('%s%s_%s_mua', obj.confSave, 'ibpc', obj.datasetLabel), Mu_a);
            dlmwrite(sprintf('%s%s_%s_taua', obj.confSave, 'ibpc', obj.datasetLabel), Tau_a);
        end
        
        function Lambda = expectedLambda(obj, Mu_a, Tau_a, Gamma)
            Lambda = ((obj.Tau_lambda.*obj.Mu_lambda) + (Gamma.*Tau_a.*Mu_a)) ...
                ./ (obj.Tau_lambda + Tau_a.*Gamma);
        end
        
        function Gamma = expectedGamma(obj, Mu_a, Tau_a, Lambda)
            Gamma = (1.5).*(Tau_a.*obj.Theta_gamma + ((Mu_a-Lambda).^2)/2)./ Tau_a;
        end
        
        function Alpha = expectedAlpha(obj, Beta, Tau_a)
            Alpha = (1 ./ obj.Theta_alpha) .* (1 + Tau_a.*Beta);
        end
        
        function Beta = expectedBeta(obj, Alpha, Tau_a)
            Beta = (1 + Alpha) .* (obj.Theta_beta + Tau_a);
        end
        
        function [Mu_a, Tau_a] = expectedPi(obj, A, T, Lambda, Gamma, Alpha, Beta)
            
            Mu_a = zeros(obj.nClasses, obj.nAgents);
            Tau_a = zeros(obj.nClasses, obj.nAgents);
            
            for j=1:obj.nClasses
            
                Tk = ones(size(A, 1), 1) *T(j, :);
                n = sum(T(j,:), 2);
                
                L = Lambda(j, :)';
                G = Gamma(j, :)';
                Al = Alpha(j, :)';
                B = Beta(j, :)';
                
                Mu_a(j,:) = (L.*G + sum(A.*Tk, 2)) ./ (G + n);
                
                AlPost = Al + n/2;
                A_hat = (sum(A.*Tk, 2) ./ n);
                if n > 0
                    BPost = B + sum((A.*Tk - (A_hat*ones(1,size(A,2)))).^2,2) ...
                        + n*G.*((A_hat-L).^2)./(2*(n+G));
                else
                    BPost = B;
                end
                Tau_a(j,:) = AlPost ./ BPost;
            end
        end
                
        function EP = expectedP(obj, T)
            EP = zeros(obj.nClasses, 1);
            for j=1:obj.nClasses
                EP(j) = obj.nu(j) + sum(T(j, :));
            end
            
            EP = EP ./ sum(EP);
        end
        
        function [ET, pT] = expectedT(obj, Tknowns, A, P, Mu_a, Tau_a)
            logPT = zeros(obj.nClasses, size(A, 2));
            pT = zeros(obj.nClasses, size(A, 2));
                                    
            for j=1:obj.nClasses                
                Mu_aj = Mu_a(j, :)' * ones(1, size(A, 2));
                Tau_aj = Tau_a(j, :)' * ones(1, size(A, 2));
                Sigma_aj = Tau_aj.^-0.5;
                
                pAiTi = normpdf(A, Mu_aj, Sigma_aj);
                                
                logPT(j, :) = log(P(j)) + sum(log(pAiTi), 1);
                pT(j, :) = exp(logPT(j,:));
            end
            pT = pT ./ (ones(obj.nClasses, 1) * sum(pT,1));            
            
            ET = pT;
            ET(:, Tknowns~=0) = 0;
            ET(sub2ind(size(logPT), Tknowns(Tknowns~=0), find(Tknowns~=0))) = 1;
        end
        
        function [ET, pT] = altPt(obj, Tknowns, C, P, Pi)
            logPT = zeros(obj.nClasses, size(C, 2));
            pT = zeros(obj.nClasses, size(C, 2));
            
            [I, J, nonZeroC] = find(C);
            nonZeroIdx = sub2ind(size(C), I, J);
            
            for j=1:obj.nClasses
                ind = sub2ind(size(Pi), ...
                    ones(size(C, 1), size(C, 2)).*j, ...
                    C, ...
                    (1:size(C, 1))'*ones(1, size(C, 2)) );
                
                logpCiTi = log( Pi(ind) );
                                                
                logPT(j, :) = log(P(j)) + sum(logpCiTi, 1);
                pT(j, :) = exp(logPT(j,:));
            end
            pT = pT ./ (ones(obj.nClasses, 1) * sum(pT,1));            
            
            ET = pT;
            ET(:, Tknowns~=0) = 0;
            ET(sub2ind(size(logPT), Tknowns(Tknowns~=0), find(Tknowns~=0))) = 1;            
        end
                
        function combinedPost = combineDecisions(obj, baseOutputs, clusters)         
            combinedPost = zeros(obj.K, size(baseOutputs, 2));
                        
            %get the relevant inputs
            if obj.K == 1
                combinedPost(1, :) = ...
                    obj.combineCluster( baseOutputs, ...
                        obj.targets);    
            else
                display('*** ignoring cluster of stupid agents: starting at k=2');
                for k=2:obj.K
                    inputPosts = baseOutputs;
                    inputPosts(clusters~=k) = 0.5;
                                        
                    combinedPost(k, :) = ...
                        obj.combineCluster( inputPosts, ...
                            obj.targets, (clusters==k));
                end
            end
            
        end
    end
    
end

