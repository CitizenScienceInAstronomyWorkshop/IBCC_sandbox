display('NEED TO ADD ERROR BARS TO SHOW DEVIATION AMONG DIFFERENT REPEATS');

nClasses = 3;
nRepeats = 20;
selectedRepeats = 1:nRepeats;%[3 9 13 14 17]; %very good outcomes for HF. Separate these out, and show against those with crappy outcomes. 

fHandles = {};
for j=1:nClasses + 4
    fHandles{j} = figure;
end

startIteration = 50;
maxNoIterations = 150;

methods = {'HFFullyDyn', 'HFStatic', 'AS', 'Random', 'OS'};%, 'HFFullyDynAll', 'HFStaticAll', 'OS'};%,  'ASAll', 'HFAS',};

for m=1:length(methods)
    selectMethod = methods{m};
    analysisDataFile = ['/homes/49/edwin/results/hcomp_sf8/' selectMethod];
%     analysisDataFile = ['/homes/49/edwin/matlab/data/hcomp_hal_arcus/analysis/hcomp2/' selectMethod];
    
    aucs_set = cell(1,nRepeats);
    Hset = cell(1,nRepeats);
    aucs_set2 = cell(1, nRepeats);
    H_set2 = cell(1, nRepeats);
    
    load(analysisDataFile); %Loads aucs_set, H_set, aucs_set2, H_set2. 
    %first attempt loads only aucs_set and H_set2 due to bug
%     aucs_set2 = aucs_set;
%     H_set = H_set2;
%     aucs_set = aucs_set(1:5);

    if exist('nRepeats','var') && length(aucs_set) > nRepeats
        aucs_set = aucs_set(1:nRepeats);
    end
    nRepeats = length(aucs_set);
    
    best = 0;
    bestIdx = 0;
    worst = nClasses;
    worstIdx = 0;
    
    longest = 0;
    for r=1:nRepeats
        
        if isempty(aucs_set{r})
            display(['Empty ' num2str(r)]);
        elseif sum(aucs_set{r}(2:end,end)) >best 
            best = sum(aucs_set{r}(2:end,end));
            bestIdx = r;
        end
        if isempty(aucs_set{r})
        elseif sum(aucs_set{r}(2:end,end))<worst
            worst = sum(aucs_set{r}(2:end,end));
            worstIdx = r;
        end            
        
        if length(aucs_set{r}) > longest
            longest = length(aucs_set{r});
        end
        if length(aucs_set{r})>maxNoIterations
            aucs_set{r} = aucs_set{r}(:,1:maxNoIterations);
            longest = maxNoIterations;
        end
    end
    display(num2str(worst));
    display(['Longest auc plot has nIterations=' num2str(longest)]);
    
    for r=1:nRepeats
        currentLength = size(aucs_set{r},2);
        if longest - currentLength > 0
            aucs_set{r} = [aucs_set{r} zeros(size(aucs_set{r},1), longest-currentLength)];
        end    
    end
    
%     nRepeats = length(aucs_set2);
%     
%     longest = 0;
%     for r=1:nRepeats
%         if length(aucs_set2{r}) > longest
%             longest = length(aucs_set2{r});
%         end
%     end
%     
%     for r=1:nRepeats
%         currentLength = size(aucs_set2{r},2);
%         if longest - currentLength > 0
%             aucs_set2{r} = [aucs_set2{r} zeros(size(aucs_set2{r},1), longest-currentLength)];
%         end    
%     end
    
    nRepeats = length(H_set);
   
    longest = 0;
    for r=1:nRepeats
        if length(H_set{r}) > longest
            longest = length(H_set{r});
        end
        if length(H_set{r})>maxNoIterations
            H_set{r} = H_set{r}(1:maxNoIterations);
            longest = maxNoIterations;
        end
    end
    
    for r=1:nRepeats
        currentLength = size(H_set{r},2);
        if longest - currentLength > 0
            H_set{r} = [H_set{r} zeros(size(H_set{r},1), longest-currentLength)];
        end    
    end
 
%     nRepeats = length(H_set2);    
%     
%     longest = 0;
%     for r=1:nRepeats
%         if length(H_set2{r}) > longest
%             longest = length(H_set2{r});
%         end
%     end
%     
%     for r=1:nRepeats
%         currentLength = size(H_set2{r},2);
%         if longest - currentLength > 0
%             H_set2{r} = [H_set2{r} zeros(size(H_set2{r},1), longest-currentLength)];
%         end    
%     end   

    nRepeats = length(aucs_set);
    
    for topic=2:size(aucs_set{1},1)
        for r=1:nRepeats

            aucs = aucs_set{r};
            if ~isempty(aucs)
                lastIdx = find(aucs(topic,:)==0);
                if ~isempty(lastIdx)
                    aucs(topic,lastIdx) = aucs(topic,lastIdx(1)-1);   
                    
                    aucs_set{r} = aucs;
                end
            end
            
%             if r <= length(aucs_set2)
%                 aucs = aucs_set2{r};
%                 if ~isempty(aucs)
%                     lastIdx = find(aucs(topic,:)==0);
%                     if ~isempty(lastIdx)
%                        aucs(topic,lastIdx) = aucs(topic,lastIdx(1)-1); 
%                         aucs_set2{r} = aucs; 
%                     end
%                 end
%             end

            if r <= length(H_set)
                H = H_set{r};

                lastIdx = find(H==0);
                if ~isempty(lastIdx) && ~isempty(H) && lastIdx(1)>1
                    H(lastIdx) = H(lastIdx(1)-1);    
                    H_set{r} = H;
                end
            end
            
%             if r <= length(H_set2)
%                 H = H_set2{r};
% 
%                 lastIdx = find(H==0);
%                 if ~isempty(lastIdx) && ~isempty(H) && lastIdx(1)>1
%                     H(lastIdx) = H(lastIdx(1)-1); 
%                     H_set2{r} = H; 
%                 end
%             end
        end
    end    
    
    sumAuc = zeros(nClasses+1, size(aucs_set{1},2));
%     sumAuc2 = zeros(nClasses+1, size(aucs_set2{1},2));
    sumH = zeros(1,size(H_set{1},2));
%     sumH2 = zeros(1,size(H_set2{1},2));
    
    divAuc = 0;
%     divAuc2 = 0;
    divH = 0;
%     divH2 = 0;
    if ~exist('selectedRepeats','var')
        selectedRepeats = 1:nRepeats;
    end
    for r=selectedRepeats
        if ~isempty(aucs_set{r})
            sumAuc = sumAuc + aucs_set{r};
            divAuc = divAuc + 1;
        end
%         if ~isempty(aucs_set2{r})
%             sumAuc2 = sumAuc2 + aucs_set2{r};
%             divAuc2 = divAuc2 + 1;
%         end
        if ~isempty(H_set{r})
            sumH = sumH + H_set{r};
            divH = divH + 1;
        end
%         if ~isempty(H_set2{r})
%             sumH2 = sumH2 + H_set2{r};    
%             divH2 = divH2 + 1;
%         end
    end
    for cl=2:nClasses+1

        %set handle 1
        figure(fHandles{cl-1});
    
%     startPoints = repmat(sumAuc(:,50), 1, size(sumAuc,2));
%     sumAuc = sumAuc - startPoints; 
%     sumAuc = sumAuc(:,50:end);
%     
%     startPoints = repmat(sumH(:,50), 1, size(sumH,2));
%     sumH = sumH - startPoints; 
%     sumH = sumH(:,50:end);    
     
        plot(sumAuc(cl,:)'./divAuc);
        hold all;
        legend(methods);
    
        title(['AUC class ' num2str(cl)])
    end
    %set handle 2 etc...
%     figure(fHandles{2});
% 
%     plot(sumAuc2(2,:)'./divAuc2);
%     hold all; 
%     legend(methods);
%     title('AUC class B')
%     
%     figure(fHandles{2});
% 
%     plot(sumAuc(3,:)'./divAuc);
%     hold all;
%     legend(methods);
%     title('AUC class B')
%     
%     figure(fHandles{4});
% 
%     plot(sumAuc2(3,:)'./divAuc2);
%     hold all; 
%     legend(methods);
%     title('AUC class D')
%     
%     figure(fHandles{3});
% 
%     plot(sumAuc(4,:)'./divAuc);
%     hold all;
%     legend(methods);
%     title('AUC class C')
%     
%     figure(fHandles{6});
% 
%     plot(sumAuc2(4,:)'./divAuc2);
%     hold all;
%     legend(methods);
%     title('AUC class F');
    
    figure(fHandles{end})
    plot(sum(aucs_set{bestIdx}(2:end,startIteration:maxNoIterations),1)./nClasses);
    hold all
    title('bestest');
    
    figure(fHandles{end-1})
    plot(sum(aucs_set{worstIdx}(2:end,startIteration:maxNoIterations),1)./nClasses);
    hold all
    title(['worstest ' selectMethod]);

    figure(fHandles{nClasses+1});
    dataToPlot = sum(sumAuc(:,startIteration:end),1)./(divAuc*nClasses);
    dataToPlot = smooth(dataToPlot);
    plot(dataToPlot);% + 0.62);
    hold all
    legend(methods);
    title('average AUC for all classes');
    
%     figure(fHandles{10});
%     
%     plot(sum(sumAuc2,1)./(divAuc2*nClasses));
%     hold all
%     legend(methods);
%     title('average AUC for all classes, second test');
    
    figure(fHandles{nClasses+2});
   
    plot(sumH(:,startIteration:end)'./divH);
    hold all;
    legend(methods);
    title('Entropy of Target Labels');  

%     figure(fHandles{6});
%     
%     plot(sumH2'./divH2);
%     hold all; 
%     legend(methods);
%     title('Entropy of Target Labels, second test');  
    
end
