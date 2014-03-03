classdef SimAgent < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        pCorrect;
        degredation = 0.002;
        id = 0;
        fired = false;
        nStable = 0;
        nResp = 0;
                
    end
    
    methods
        function obj = SimAgent(id, degrade, nIterations, pCorrect, nClasses, nStable)
            
            if ~exist('nClasses','var')
                display('assuming 2 class problem for sim agent');
                nClasses = 2;
            end
            
            if exist('pCorrect','var')
                obj.pCorrect = pCorrect;
            else
                obj.pCorrect = (rand(1,1) .* 0.4) + 1/nClasses;
            end
            
            if exist('nStable', 'var')
                obj.nStable = nStable;
            end
            
            if exist('nIterations','var')    
                if nIterations > 0
                    obj.degredation = (obj.pCorrect - (1/nClasses))/nIterations;
                else
                    display('degrade to always clicking 1');
                    obj.degredation = obj.pCorrect + 1;
                end
            end
            if degrade==0
                obj.degredation = 0;
            elseif degrade==-1
                obj.degredation = -obj.degredation;
            end
            obj.id = id;
        end
        
        function answer = getResponse(obj, truth, doc)
            
            display(['worker ' num2str(obj.id) ' has pcorrect ' num2str(obj.pCorrect)]);
            
            if obj.pCorrect<=0
                display('I am just going to lazily click category 1 without reading the document...');
                answer = 1;
                return
            end
            
            ti = truth(doc,:);
            ti = find(ti);
            if isempty(ti)
                ti = 1;
            else
%                 display('true positive');
            end
            ti = ti(randi(length(ti),1,1)); %where multiple classes are present, pick at random
            
            r = rand(1,1);
            if r>obj.pCorrect
                %gets it wrong, usually defaulting to "none-of-the-above" class 1
%                 if ti~=1
%                     ti = 1;
%                 else
%                     ti = 1 + randi(size(truth,2)-1, 1,1);
%                 end

                noTopic = rand(1,1)>0.7;
                if noTopic %does the agent miss the topic of this document? If it really is none of the above, the label is changed randomly
                    answer = 1;
                else
                    answer = ti;
                end
                
                while answer==ti
                    answer = randi(size(truth,2),1,1);
                end
            else
                answer = ti;
            end
            
            minPCorrect = 1/size(truth,2);
            
            if obj.pCorrect > minPCorrect && obj.nResp > obj.nStable
                obj.pCorrect = obj.pCorrect - obj.degredation;
                if obj.pCorrect < minPCorrect %catch it when it dips below 0.5
                    obj.pCorrect = 0;%minPCorrect;
                end
            elseif obj.nResp > obj.nStable
                display(['agent ' num2str(obj.id) ' is stupid and should be fired']);
            end
            
            obj.nResp = obj.nResp + 1;
            
            display(['truth: ' num2str(ti) ', answer: ' num2str(answer)]);
        end
    end
    
end

