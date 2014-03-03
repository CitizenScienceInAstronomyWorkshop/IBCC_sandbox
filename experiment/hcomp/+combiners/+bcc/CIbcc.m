classdef  CIbcc < combiners.bcc.Ibcc
 %Ibcc for continuous label values. In future should handle a mixture of
 %continuous and discreet labels correctly. Discreet labels should also
 %include "don't know" or correspond to probabilities of each class (e.g
 %rate between 1 and 5 your belief in class 1).
    properties
    end
    
    methods (Static)
        function sl = shortLabel
            sl = 'CIbcc';
        end
           
    end
    
    methods       
        function obj = CIbcc(nAgents, K, targets, agents, nClasses)
            
            if nargin < 5
                nClasses = 2;
            end 
            
            obj@combiners.bcc.Ibcc(nAgents, K, targets, agents, nClasses);
      
            obj.label = 'Continuous Indep. Bayesian Classifier Combination';
            
        end
         
%         function BetaSample = sampleBeta(obj, Alpha)     
%             %Beta is the set of parameters for the confusion distribution.
%             %Confusion is no longer a matrix of fixed values as the
%             %classifier outputs are continuous between 0 and 1.
%         end
        
%         %Find the p(Ti | Ci, P, Conf, Alpha) =  p(Ti | Ci, P, Conf)
%         function postTi = posteriorTi(obj, C, Conf, P, members)
% 
%         end
        
        function C = prepareC(obj, post)
            C = post;
        end
        
        function combinedPost = combineCluster(obj, post, T, members)          
           
            Tknown = T(1:obj.nKnown);
            C = obj.prepareC(post);        
            obj.nClassifiers = size(C,1);
            
            if nargin > 3
                gibbsSampler = combiners.bcc.IbccSampling.GibbsSampling(obj, members);
            else
                gibbsSampler = combiners.bcc.IbccSampling.GibbsSampling(obj);
            end
            
            expectedT = gibbsSampler.sample(C, Tknown);
%             
%             if nargin > 3
%                 expectedT = obj.sample(C, Tknown, members);
%             else
%                 expectedT = obj.sample(C, Tknown);
%             end
            combinedPost = expectedT;
        end  
        
        function expectedPost = sample(obj, C, Tknown)
                        
            nSamples = size(C, 2);
                        
            [Alpha, P, T, Tknown] = obj.initVariables(Tknown, nSamples);
            
            currentInterval = obj.sampleInterval;
            nSamplesCollected = 0;
            
            expectedPost = zeros(1, nSamples);
            
            i = 0;            
            while nSamplesCollected < obj.gibbsSamples
                i = i+1;
                display(['gibbs samples collected: ' num2str(nSamplesCollected) ', ' num2str(i)]);

                %Conf conditional is prop. to product of dirichlets
                %p(pi|alpha)p(pi|c,t) = p(pi| alpha + confCounts)
                confCounts = obj.confusionCounts(T, C, nSamples);
                Conf = obj.combiner.sampleConf(Alpha + confCounts);
                if Conf(1, 2, 1) > 0.8 || Conf(1, 2, 2) > 0.8 
                    display('wrong conf');
                end
                %uses lambda as prior
                Alpha = obj.sampleAlpha(Conf, Alpha); 
                
                %need to adapt this so that we include some correct values
                %of T. What about combination on the fly?
                [Tsample, T] = obj.combiner.sampleT(C, Conf, P, Tknown, obj.members);
                
                P = obj.sampleP(T); %also uses nu as a prior
                                
                if i>obj.burnIn && currentInterval >= obj.sampleInterval
                    %T records labels as 1 or 2
                    expectedPost = expectedPost + T - 1;
                    currentInterval = 1;
                    nSamplesCollected = nSamplesCollected + 1;
                else
                    currentInterval = currentInterval + 1;                    
                end
            end  
            expectedPost = expectedPost ./ nSamplesCollected;
        end        
    end
    
end

