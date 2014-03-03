classdef  IbpcVb2 < combiners.bcc.Ibcc
    %Independent Bayesian Probability Combination: combines
    %probabilisitic outputs from classifiers (rather than discrete
    %classifications

    properties
        maxIt = 500;
        convIt = 10;
                
        %Extra Hyperhyperparameters
        A
        B
    end
 
    methods (Static)
        function sl = shortLabel
            sl = 'Ibpc2';
        end   
    end
    
    methods       
        function obj = IbpcVb2(nAgents, K, targets, agents, nClasses, nScores, ...
                lambdaMag, lambdaSym, nu)
            %hyperparameter priors
            if nargin < 7
                lambdaMag = 1;
                lambdaSym = 0.2;
            end

            obj@combiners.bcc.Ibcc(nAgents, K, targets, agents, nClasses, nScores, ...
                lambdaMag, lambdaSym, nu);

            obj.label = 'I.B.P.C. 2'; 

            obj.nonDet = false;
            
            if nClasses>2
                display('Not suitable for multiclass yet');
            end
                        
            offDiag = lambdaMag * (1-lambdaSym);
            diag = lambdaMag * lambdaSym;
                        
            obj.A = ones(nClasses, nScores, nAgents) * offDiag;
            obj.A(sub2ind(size(obj.A), 1:obj.nClasses, 1:obj.nClasses)) ...
                    = ones(obj.nClasses, 1) * diag;

            obj.B = ones(nClasses, nScores, nAgents) * offDiag;   
            obj.B(sub2ind(size(obj.B), 1:obj.nClasses, 1:obj.nClasses)) ...
                    = ones(obj.nClasses, 1) * diag;
            
            for k=1:obj.nClassifiers
                obj.A(:, :, k) = obj.A(:, :, 1);
                obj.B(:, :, k) = obj.B(:, :, 1);
            end
        end
        
        function Post = preparePost(obj, BasePost)
            nSamples = size(BasePost, 2);
            
            Post = zeros(obj.nScores, nSamples, obj.nClassifiers);
            %not suitable for >2 classes as posteriors provided by base
            %classifiers are just a single value
            Post(1, :, :) = reshape(BasePost', 1, nSamples, obj.nClassifiers);
            Post(2, :, :) = 1 - Post(1, :, :);
        end
        
        function [Lambda, Alpha, P, T, Tknown] = initVariables(obj, Tknown, Post)
                     
            Lambda = obj.A .* obj.B;
            
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
            
            Alpha = obj.expectedAlpha(Post(:, Tknown>0), T(:, Tknown>0), Lambda);
        end        
                
        function combinedPost = combineCluster(obj, Post, T, members)          
            obj.nKnown = round(obj.nKnown);
            Tknown = T;
            
            obj.nClassifiers = size(Post,1);
            Post = obj.preparePost(Post);
            
            converged = false;
            pT = 0;
            
            nIt = 0;
            cIt = 0;
            
            %change T into the matrix format for multi-class problems
            [Lambda, Alpha, P, T, Tknowns] = ...
                obj.initVariables(Tknown, Post);
            while ~converged
                oldPT = pT;
                
                [T, pT] = obj.expectedT(Tknowns, Post, P, Alpha);
                
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
                    display(sprintf('IbPcVB2 iterations: %d, iterations while converged %d, ln(pt)=%f', nIt, cIt, log(pT)));
                end
                
                %calc expected values
                Lambda = obj.expectedLambda(Alpha, size(Post, 2));
                Alpha = obj.expectedAlpha(Post, T, Lambda);
                P = obj.expectedP(T);
                                
                if pT < oldPT 
                    display(sprintf('IbpcVb2: %d, %f', nIt, log(pT)));
                end
            end
            
            %need to change this for >2 classes as it only gives
            %probability of T==2
            combinedPost = T(2, :);
            
            %save the expected confusion distributions
            dlmwrite(sprintf('%s%s_%s_alpha', obj.confSave, 'ibpc2', obj.datasetLabel), Alpha);
        end
        
        function Lambda = expectedLambda(obj, Alpha, n)
            postA = obj.A + n;
            postB = obj.B + n.*Alpha;
            Lambda = postA .* postB;
        end
                 
        function Alpha = expectedAlpha(obj, Post, T, Lambda)
            Alpha = zeros(obj.nClasses, obj.nScores, obj.nClassifiers);
            for j=1:obj.nClasses
                for p=1:obj.nScores
                    for k=1:obj.nClassifiers
                        if size(Post, 2) == 0
                            sumLogData = 0;
                        else
                            Data = Post(p, :, k) .* T(j, :);
                            sumLogData = sum(log(Data), 2); %
                            display('might have accidentally deleted some code here - cannot remember what');
                        end
                        Alpha(j, p, k) = sumLogData + (1) ./ Lambda(j, p, k);
                    end
                end
            end
        end
                
        function EP = expectedP(obj, T)
            EP = zeros(obj.nClasses, 1);
            for j=1:obj.nClasses
                EP(j) = obj.nu(j) + sum(T(j, :));
            end
            
            EP = EP ./ sum(EP);
        end
        function f = logDirPdf(obj, X, Alpha)
            logBeta = sum(log(gamma(Alpha)), 2) - log(gamma(sum(Alpha, 2)));
            
            f = sum(log(X).*(Alpha-1), 2) - logBeta;
        end
        
        function [ET, pT] = expectedT(obj, Tknowns, Post, P, Alpha)
            logPT = zeros(obj.nClasses, size(Post, 2));
            pT = zeros(obj.nClasses, size(Post, 2));
                                    
            for j=1:obj.nClasses                
                Alpha_j = Alpha(j, :, :);
                
                %Inputs to dirpdf: set of values, obj.nScores
                
                Alpha_j = reshape(Alpha_j, obj.nScores, obj.nClassifiers)';
                
                logPPostiTi = zeros(obj.nClassifiers, size(Post, 2));
                
                for n=1:size(Post, 2)
                    Post_n = reshape(Post(:, n, :), size(Post, 1), size(Post, 3))';
                    logPPostiTi(:, n) = obj.logDirPdf(Post_n, Alpha_j);
                end
                                
                logPT(j, :) = log(P(j)) + sum(logPPostiTi, 1);
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

