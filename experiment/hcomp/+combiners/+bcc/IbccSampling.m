classdef IbccSampling < combiners.bcc.Ibcc
    %IBCCSAMPLING IBCC class that uses sampling
    %   Detailed explanation goes here
    
    properties
        lambda = [0.5 1; 1 0.5];%[0.5 1; 1 0.5];%[0.4 1; 1 0.4]; %[1 1; 1 1];
        lambdaMag; lambdaSym;        
    end
    
    methods (Static)
        function sl = shortLabel
            sl = 'IbccSampling';
        end
    end
     
    methods
        
        function obj = IbccSampling(bccSettings, nAgents, K, targets, agents, nClasses, nScores, nu)
            obj@combiners.bcc.Ibcc(bccSettings, nAgents, K, targets, agents, nClasses, nScores);
            
            while size(obj.lambda, 2) < nScores
                extension = [0.5; 0.5];
                obj.lambda = [obj.lambda extension]; 
            end
            
            %hyperparameter priors            
            obj.lambdaMag = bccSettings.lambdaMag(bccSettings.iLambdaMag); 
            obj.lambdaSym = bccSettings.lambdaSym(bccSettings.iLambdaSym);                            
            diagLambda = bccSettings.lambdaMag(1) * bccSettings.lambdaSym(1);
            offDiagLambda = bccSettings.lambdaMag(1) * (1-bccSettings.lambdaSym(1));
            
            obj.lambda = offDiagLambda .* ones(obj.nClasses, obj.nScores);
            obj.lambda(sub2ind(size(obj.lambda), 1:obj.nClasses, 1:obj.nClasses)) ...
                = diagLambda .* ones(obj.nClasses, 1);
            
            obj.lambda = reshape(obj.lambda, size(obj.lambda, 1)*size(obj.lambda, 2), 1);
            obj.lambda = obj.lambda * ones(1, obj.nAgents);
            obj.lambda = reshape(obj.lambda, obj.nClasses, obj.nScores, obj.nAgents);
            
            if nargin > 7
                obj.Nu0 = nu;
            end
            
            obj.setTrustedAgents();
                        
            obj.combinerInfo = ['lambda(:,:,1) ' mat2str(obj.lambda(:,:,1)) ', nu ' mat2str(obj.Nu0)];            
            obj.nonDet = true; %non-deterministic
        end        
        
        function setTrustedAgents(obj, trustFinalAgent)
            if nargin > 1
                obj.trustFinalAgent = trustFinalAgent;
            end
            for trusted=obj.trustFinalAgent
                trustedStart = obj.nAgents - length(obj.trustFinalAgent);
                if isempty(obj.trustedLambda)
                    obj.lambda(:,:,trustedStart+trusted) = obj.lambdaMag(1+trusted)*obj.lambdaSym(1+trusted);
                    idxs = sub2ind(size(obj.lambda), ...
                        1:obj.nClasses, 1:obj.nClasses, ones(1,obj.nClasses)*trustedStart+trusted);
                    obj.lambda(idxs) = obj.lambdaMag(1+trusted)*(1-obj.lambdaSym(1+trusted));
                else
                    obj.lambda(:, :, trustedStart+trusted) = obj.trustedLambda(:,:,trusted);
                end
            end
        end   
        
                
        function combinedPost = combineCluster(obj, post, T, members)          
           
            obj.nKnown = round(obj.nKnown);
            Tknown = T;
            
            C = obj.prepareC(post);
                          
            %gibbs generates a bunch of samples, and takes the mean as
            %the expected value of the distribution.
            if strcmp(obj.samplingType, 'Gibbs')
                if nargin > 3
                    gibbsSampler = combiners.bcc.ibccsampling.GibbsSampling(obj, obj.nAgents, members);
                else
                    gibbsSampler = combiners.bcc.ibccsampling.GibbsSampling(obj, obj.nAgents);
                end
                expectedT = nan;
                i = 1;
                while isnan(expectedT)
                    if i > obj.maxSampleAttempts
                        break
                    end
                    [expectedT, expectedConf, expectedAlpha] = gibbsSampler.sample(C, Tknown);
                    i = i+1;
                end    
                %save the expected confusion matrix
                %dlmwrite(sprintf('%s%s_%s', obj.confSave, 'ibccgibbs',
                %obj.datasetLabel), expectedConf);
                combinedPost = expectedT;
            elseif strcmp(obj.samplingType, 'MH')
                if nargin > 3
                    mhSampler = combiners.bcc.ibccsampling.MHSampling(obj, members);
                else
                    mhSampler = combiners.bcc.ibccsampling.MHSampling(obj);
                end
                combinedPost = mhSampler.sample(C, Tknown);
            end
            
            obj.Alpha = obj.Alpha + expectedAlpha;
        end
               
    end
    
end

