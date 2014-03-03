if ~exist('keepAgents','var') || keepAgents == false || ~exist('agents','var')
    agents = cell(nWorkers, 1);

%     agents{1} = SimAgent(1, false, 0, 1/nClasses, nClasses);
%     agents{2} = SimAgent(2, false, 0, 1/nClasses, nClasses);
%     agents{3} = SimAgent(3, false, 0, 1/nClasses, nClasses);
%     agents{4} = SimAgent(4, false, 0, 1/nClasses, nClasses);
%     agents{5} = SimAgent(5, false, 0, 1/nClasses, nClasses);
%
%     agents{1} = SimAgent(1, true, nToCollect/(nWorkers*4), 0.99, nClasses, 10);
%     agents{2} = SimAgent(2, true, nToCollect/(nWorkers*4), 0.5, nClasses, 15);
%     agents{4} = SimAgent(4, true, nToCollect/(nWorkers*2), 0.8, nClasses, 3);
%     agents{5} = SimAgent(5, true, nToCollect/(nWorkers*4), 0.5, nClasses, 10);
%     agents{3} = SimAgent(3, true, nToCollect/(nWorkers*2), 0.3, nClasses, 5);

    %USE THIS FOR REAL EXPERIMENTS
%  %      agents{1} = SimAgent(1, true, 1, 1/nClasses^0.07, nClasses, 10);%24);%2*18);
%  %      agents{2} = SimAgent(2, true, 1, 1/nClasses^0.33 , nClasses, 20);
%  %      agents{3} = SimAgent(3, true, 1, 1/nClasses, nClasses, 20);

%%%%%%%%

    agents{1} = SimAgent(1, true, 1, 0.95, nClasses, 10);%24);%2*18);
    agents{2} = SimAgent(2, true, 1, 0.8 , nClasses, 15);
    agents{3} = SimAgent(3, true, 1, 0.5, nClasses, 20);

    agents{4} = SimAgent(4, true, 1, 0.95, nClasses, 20);
    agents{5} = SimAgent(5, true, 1, 0.8 , nClasses, 25);
%     agents{6} = SimAgent(6, true, 1, 0.5, nClasses, 200);
    
%     agents{7} = SimAgent(7, true, 0, 0.95, nClasses, 300);
%     agents{8} = SimAgent(8, true, 0, 0.79 , nClasses, 350);
%     agents{9} = SimAgent(9, true, 0, 0.5, nClasses, 200);

%%%%%%%%

%      agents{10} = SimAgent(10, true, 1, 0.5, nClasses, 20);

%     agents{4} = SimAgent(4, true, 1, 0.8, nClasses, 40);
%     agents{5} = SimAgent(5, true, 1, 0.55, nClasses, 28);
    
%     %USE THIS FOR SEEING IF WE GET A GOOD AUC
%     agents{1} = SimAgent(1, true, 1, 0.99, nClasses, 2*180);
%     agents{2} = SimAgent(2, true, 1, 0.55, nClasses, 2*150);
%     agents{4} = SimAgent(4, true, 1, 0.8, nClasses, 2*200);
%     agents{5} = SimAgent(5, true, 1, 0.55, nClasses, 2*130);
%     agents{3} = SimAgent(3, true, 1, 0.25, nClasses, 2*100);
    

%     agents{1} = SimAgent(1, true, 1, 0.99, nClasses, 180);
%     agents{2} = SimAgent(2, true, 1, 0.5, nClasses, 150);
%     agents{4} = SimAgent(4, true, 1, 0.8, nClasses, 200);
%     agents{5} = SimAgent(5, true, 1, 0.5, nClasses, 120);
%     agents{3} = SimAgent(3, true, 1, 0.3, nClasses, 100);


%
%     for k=6:nWorkers-nDegrade
%         agents{k} = SimAgent(k, false);
%     end
%     if(numel(k)==0)
%         k=6;
%     end
%     for k=(k+1):nWorkers
%         agents{k} = SimAgent(k, true, nToCollect/(nWorkers*2));
%     end

    if exist('initialLabels','var') && ~isempty(initialLabels)
        %Make sure the agents have deteriorated correctly
        getSimResponses(agents, initialLabels(:,1), ...
            initialLabels(:,2), [], qRels(:,[1; chosenIdx]) );
    end
else
    display('keeping same agents - if agents have degraded they should be reset');
end