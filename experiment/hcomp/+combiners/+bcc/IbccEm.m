classdef  IbccEm < combiners.bcc.Ibcc
 %T should be in matrix format with rows corresponding to classes, and
 %entries giving the probability of a class
    
 properties
    P
    T
    maxIt = 200;  
    
    totalIts = 0;
 end
 
    methods (Static)
        function sl = shortLabel
            sl = 'IbccEm';
        end   
    end
    
    methods       
        function obj = IbccEm(nAgents, K, targets, agents, nClasses, nScores)
            obj@combiners.bcc.Ibcc(nAgents, K, targets, agents, nClasses, nScores);

            obj.label = 'EM I.C.C.'; 
            obj.nonDet = false;
        end
        
        function setTrustedAgents(obj, trustFinalAgent)
            %do nothing as not priors in this version
        end
        
        function [Pi, P, T, Tknown] = initVariables(obj, Tknown, C)
            P = [0.5 0.5];           
            
            %Tknown = Tknown + 1;
            T = zeros( obj.nClasses, length(Tknown));
            for n=1:length(Tknown)
                %leave the unknown ones as equal probability of each
                if Tknown(n) == 0
                    T( :, n) = 1 ./ obj.nClasses;
                else
                    T( Tknown(n), n) = 1;
                end
            end
            
            Pi = obj.expectedPi( C, T);
        end        
                
        function combinedPost = combineCluster(obj, post, T, members)          
            obj.nKnown = round(obj.nKnown);
            Tknown = T;
            C = obj.prepareC(post);
            converged = false;
            %likelihood we wish to maximise
            L = 0;
            
            %change T into the matrix format for multi-class problems
            [Pi, P, T, Tknowns] = obj.initVariables(Tknown, C);
            
            nIt = 0;
            %bestT = T;
            convThreshold = 10^-8;
            ET = T;
            
            while ~converged                
                %calculate expected values of latent variables
                oldET = ET;
                ET = obj.expectedT(Tknowns, C, log(P), log(Pi));
                
                %evaluate likelihood
                oldL = L;
                L = obj.completeDataLogLikelihood(C, ET, Pi, P);
                   
                nIt = nIt + 1;
                if sum(abs(ET(2,:)-oldET(2,:))) < convThreshold || nIt > obj.maxIt
                    converged = true;
                    display(sprintf('EM iterations: %d', nIt));
                else
                    %maximise L wrt Pi and P
                    Pi = obj.expectedPi(C, ET);
                    P = obj.expectedP(ET);
                end
                
                if L <= oldL 
                    %display(sprintf('ICC EM: %d, %f', nIt, L));
                else
                    %bestT = ET;
                end
            end
            
            %ET = bestT;
            %need to change this for >2 classes as it only gives
            %probability of T==2
            combinedPost = ET(2, :);
            
            obj.P = P;
            obj.Pi = Pi;
            
            obj.T = ET;
            
            %save the expected confusion matrix
            if ~isempty(obj.confSave)
                %dlmwrite(sprintf('%s%s_%s', obj.confSave, 'ibccem', obj.datasetLabel), Pi);
            end
            
            nItFile = sprintf('%s%s', obj.confSave, 'nIts_em');
            if exist(nItFile, 'file');
                prevNIts = dlmread(nItFile);
                prevNIts = [prevNIts nIt];
            else
                prevNIts = nIt;
            end
            dlmwrite(nItFile, prevNIts);
        end
        
        function L = completeDataLogLikelihood(obj, C, ET, Pi, P)
            %L = sum(sum(ET .* log(ET), 1), 2);
            
            if obj.scoreSet
                nSamples = size(ET, 2);
                logJoint = sparse(obj.nClasses, nSamples);
                try
                    validDecs = C{3}(C{3}~=0);
                    validN = C{2}(C{3}~=0)';
                catch me
                    display(['index for E[T] failed: ' me.message]);
                end
                for j=1:obj.nClasses
                
                    indx = sub2ind(size(Pi), ...
                                j*ones(numel(validDecs), 1), ...
                                validDecs, ...
                                C{1}(C{3}~=0));

                    lnPCT = sparse(j*ones(1, length(validN)), validN, ...
                        log(Pi(indx)), obj.nClasses, nSamples);
                    try
                        logJoint = logJoint + lnPCT;
                        logJoint(j,:) = logJoint(j,:) + log(P(j));
                    catch me
                        display('log joint failed in ICC EM');
                    end
                end
                
                L = sum(sum(ET .* logJoint, 2), 1);
            else
                L = 0;

                [I, J, nonZeroC] = find(C);

                if size(I, 2) > size(I, 1)
                    I = I';
                end
                if size(nonZeroC, 2) > size(nonZeroC, 1)
                    nonZeroC = nonZeroC';
                end

                for j=1:obj.nClasses
                    confIdxs = sub2ind(size(Pi), ones(length(I), 1).*j, nonZeroC, I);

                    logConfs = log(Pi(confIdxs));
                    sumLogConfs = sum(logConfs, 1);
                    logJointSampleProb = log(P(j)) + sumLogConfs; %this used to be multiplied and worked but seems wrong.

                    L = L + sum(ET(j, :) .* logJointSampleProb, 2);
                end
            end
        end
  
        function Pi = expectedPi(obj, C, T)
            Pi = zeros(obj.nClasses, obj.nScores, obj.nAgents);
            
            %dupMatrix = ones(obj.nAgents, 1);
            
            Count = obj.voteCounts(C,T);
            
            for j=1:obj.nClasses
                
                %Tj = dupMatrix * T(j, :);
                
                for l=1:obj.nScores 
                    Countjl = Count(j,l,:);%sum( (C==l) .* ( Tj ), 2);
                    Pi(j, l, :) = Countjl;
                end
                
                PiSum = sum(Pi(j, :, :), 2);
                
                for l=1:obj.nScores
                    Pi(j, l, :) = Pi(j, l, :) ./ PiSum;
                    Pi(j, l, PiSum==0) = 1 ./ obj.nScores;
                end
            end
        end
                
        function EP = expectedP(obj, T)
            EP = zeros(obj.nClasses, 1);
            for j=1:obj.nClasses
                EP(j) = sum(T(j, :));
            end
            
            EP = EP ./ sum(EP);
        end

    end
    
end

