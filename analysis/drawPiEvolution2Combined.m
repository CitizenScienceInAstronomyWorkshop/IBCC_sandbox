function [figHandleInc figHandleDec figHandleIncdec] = drawPiEvolution2Combined(Alpha, Kt, changePoints, figHandleInc, figHandleDec, figHandleIncdec, trans)
%Kt is the list of classifier IDs for each sample in the sample order, i.e. C{1} from the scoreset passed to the combiners
%drawPiEvolution Draw a 3-D graph showing the evolution of Pi over time.
%Works only when there are 3 possible scores (e.g. GZSN). A separate
%graph is drawn for each target label value.

if ~exist('trans','var')
    trans = false;
end

minN = 0;

Kunique = unique(Kt); %the ids of the base classifiers


normTerm = repmat(sum(Alpha,2), [1 size(Alpha,2) 1]);
Pi = Alpha ./ normTerm;


if size(Pi, 2) ~=2
    display('Pi Evolution 2scores is only implemented for confusion matrices with 2 possible scores.');
end

% Kfiltered = [];

Kinc = [];
Kdec = [];
Kincdec = [];

for k=Kunique'
    
    Nk = find(Kt==k); %the data points corresponding to this base classifier
    
    nNk = length(Nk);
    if nNk < minN
        continue;
    end
    if k>length(changePoints)
        continue; 
        %this agent didn't change so we won't bother drawing it
    elseif ismember(k, [])%[2 3 4 5 6])
        %decreasing because we always made the odd-numbered agents decrease
        %in the simulation
        Kdec = [Kdec k];
    elseif ismember(k,[1 2 13 14 15])
        Kinc = [Kinc k];
    elseif ismember(k, [1])% 3 5 7 8 9 10 11 12])
        Kincdec = [Kincdec k];
    end
end

display(num2str(length(Kinc)));
display(num2str(length(Kdec)));

if ~exist('figHandleInc','var') || isempty(figHandleInc)
    figHandleInc=figure('Position', [400 400 1000 700]); 
else
    set(0, 'currentfigure', figHandleInc);
end
%%% 1. SEPARATE AGENTS IMPROVING AND DECREASING
%%% 2. CENTRE THE X-DATA on the changepoints
    colours = {[0,0,0], [1,0,0]};
    maxX = 0;
    minX = 0;
    for j=1:size(Pi,1)
        plot([-1 -1],'Color', colours{j});hold all
    end
    for k=Kinc
        Nk = find(Kt==k); %the data points corresponding to this base classifier
        
        display(['number of data points for ' num2str(k) ': ' num2str(length(Nk))]);
        hold all

        legendStrings = cell(size(Pi,1), 1);

        X = zeros(size(Pi,1), length(Nk));
        Z = zeros(size(Pi,1), length(Nk));    

        for j=1:size(Pi,1)

            Pij = Pi(j,:,:);%Alpha(j,:,:);%

            set(gca, 'FontSize', 18);

            for n=Nk'
                Pijt = Pij(1, :, n);
    %             dist = 1-Pijt;
                x = Pijt(2);%./(Pijt(1)+Pijt(2));
                z = find(Nk==n);

                X(j, z) = x;
                Z(j, z) = z - changePoints(k);
                if max(z - changePoints(k))>maxX
                    maxX = max(z - changePoints(k));
                end     
                if min(z - changePoints(k))<minX
                    minX = min(z- changePoints(k));
                end
            end    
            w = 50;
            if ~trans
                plot(Z(j,:), X(j,:), '-', 'linewidth', 2, 'color', colours{j});  
            else
                patchline(Z(j,:), X(j,:), 'linewidth', 2, 'EdgeColor', colours{j}, 'EdgeAlpha', 0.4);  
            end
        end

%         textYOffset = -0.05;
%         textXOffset = -0.07; 
%         h1 = text(0.5+textXOffset, textYOffset, 'score=0');
%         h2 = text(0.5+textXOffset, 1-textYOffset, 'score=1');

        for j=1:size(Pi,1)
            legendStrings{j} = ['\pi^k_' num2str(j-1)];    
        end

        legend(legendStrings);        
%         set(gca,'ZColor','w');
%         set(gca,'ZTick', []);        
%         set(gca, 'YTick', []);
    %     set(gca, 'XTickLabel', []);
%         set(h1, 'FontSize', 18);
%         set(h2, 'FontSize', 18);
        xlabel('Sample No. Relative to Change-Point');
        ylabel('Probability of Response=1');
        
        axis([minX-50 maxX+50 0 1]);

    %     saveas(h, sprintf('/homes/49/edwin/matlab/combination/results/dynamicVB/paper/single_dynamics_synth/pi2d_%d.png', k), 'png');    
    %     saveas(h, sprintf('/homes/49/edwin/matlab/combination/results/dynamicVB/paper/single_dynamics_synth/pi2d_%d.fig', k), 'fig');    
    end
    
    plot([0 0], [0 1], 'Color', 'black');
    
    title('Evolution of E[\pi_j^k] for Improving Agents');
    
    grid on
    
if ~exist('figHandleIncDec','var') || isempty(figHandleIncdec)
    figHandleIncdec=figure('Position', [400 400 1000 700]); 
else
    set(0, 'currentfigure', figHandleIncdec);
end
%%% 1. SEPARATE AGENTS IMPROVING AND DECREASING
%%% 2. CENTRE THE X-DATA on the changepoints

    maxX = 0;
    minX = 0;
    for j=1:size(Pi,1)
        plot([-1 -1],'Color', colours{j});hold all
    end
    for k=Kincdec
        Nk = find(Kt==k); %the data points corresponding to this base classifier
        
        display(['number of data points for ' num2str(k) ': ' num2str(length(Nk))]);
        hold all

        legendStrings = cell(size(Pi,1), 1);

        colours = {[0,0,0], [1,0,0]};

        X = zeros(size(Pi,1), length(Nk));
        Z = zeros(size(Pi,1), length(Nk));    

        for j=1:size(Pi,1)

            Pij = Pi(j,:,:);

            set(gca, 'FontSize', 18);

            for n=Nk'
                Pijt = Pij(1, :, n);
    %             dist = 1-Pijt;
                x = Pijt(2)./(Pijt(1)+Pijt(2));
                z = find(Nk==n);

                X(j, z) = x;
                Z(j, z) = z - changePoints(k);
                if max(z - changePoints(k))>maxX
                    maxX = max(z - changePoints(k));
                end
                if min(z - changePoints(k))<minX
                    minX = min(z - changePoints(k));
                end
            end    
            w = 50;
            if ~trans
                plot(Z(j,:), X(j,:), '-', 'linewidth', 2, 'color', colours{j});  
            else
                patchline(Z(j,:), X(j,:), 'linewidth', 2, 'EdgeColor', colours{j}, 'EdgeAlpha', 0.4);  
            end
        end

%         textYOffset = -0.05;
%         textXOffset = -0.07;
%         h1 = text(0.5+textXOffset, textYOffset, 'score=0');
%         h2 = text(0.5+textXOffset, 1-textYOffset, 'score=1');

        for j=1:size(Pi,1)
            legendStrings{j} = ['\pi^k_' num2str(j-1)];    
        end

        legend(legendStrings);        
%         set(gca,'ZColor','w');
%         set(gca,'ZTick', []);        
%         set(gca, 'YTick', []);
    %     set(gca, 'XTickLabel', []);
%         set(h1, 'FontSize', 18);
%         set(h2, 'FontSize', 18);
        xlabel('Sample No. Relative to Change-Point');
        ylabel('Probability of Response=1');
        
        axis([minX-50 maxX+50 0 1]);

    %     saveas(h, sprintf('/homes/49/edwin/matlab/combination/results/dynamicVB/paper/single_dynamics_synth/pi2d_%d.png', k), 'png');    
    %     saveas(h, sprintf('/homes/49/edwin/matlab/combination/results/dynamicVB/paper/single_dynamics_synth/pi2d_%d.fig', k), 'fig');    
    end
    
    plot([0 0], [0 1], 'Color', 'black');
    title('Evolution of E[\pi_j^k] for Agents that Improve then Deteriorate');    
    grid on
    
    if ~exist('figHandleDec','var')|| isempty(figHandleDec)
        figHandleDec=figure('Position', [400 400 1000 700]); 
    end
    set(0, 'currentfigure', figHandleDec);
    
    maxX = 0;
    minX = 0;
    for j=1:size(Pi,1)
        plot([-1 -1],'Color', colours{j});hold all
    end    
    for k=Kdec
        Nk = find(Kt==k); %the data points corresponding to this base classifier

        display(['number of data points for ' num2str(k) ': ' num2str(length(Nk))]);
        hold all

        legendStrings = cell(size(Pi,1), 1);

        colours = {[0,0,0], [1,0,0]};

        X = zeros(size(Pi,1), length(Nk));
        Z = zeros(size(Pi,1), length(Nk));    
        changePoints(k)
        for j=1:size(Pi,1)

            Pij = Pi(j,:,:);

            set(gca, 'FontSize', 18);

            for n=Nk'
                Pijt = Pij(1, :, n);
    %             dist = 1-Pijt;
                x = Pijt(2)./(Pijt(1)+Pijt(2));
                z = find(Nk==n);
                
                
                
                X(j, z) = x;
                Z(j, z) = z - changePoints(k);
                if max(z - changePoints(k))>maxX
                    maxX = max(z - changePoints(k));
                end         
                if min(z - changePoints(k))<minX
                    minX = min(z - changePoints(k));
                end                
            end    
            w = 50;
            if ~trans
                plot(Z(j,:), X(j,:), '-', 'linewidth', 2, 'color', colours{j});  
            else
                patchline(Z(j,:), X(j,:), 'linewidth', 2, 'EdgeColor', colours{j}, 'EdgeAlpha', 0.4);  
            end
        end

%         textYOffset = -0.05;
%         textXOffset = -0.07;
%         h1 = text(0.5+textXOffset, textYOffset, 'score=0');
%         h2 = text(0.5+textXOffset, 1-textYOffset, 'score=1');

        for j=1:size(Pi,1)
            legendStrings{j} = ['\pi^k_' num2str(j-1)];    
        end

        legend(legendStrings);        
%         set(gca,'ZColor','w');
%         set(gca,'ZTick', []);        
%         set(gca, 'YTick', []);
    %     set(gca, 'XTickLabel', []);
%         set(h1, 'FontSize', 18);
%         set(h2, 'FontSize', 18);
        xlabel('Sample No. Relative to Change-Point');
        ylabel('Probability of Response=1');
        
        axis([minX-50 maxX+50 0 1]);

    %     saveas(h, sprintf('/homes/49/edwin/matlab/combination/results/dynamicVB/paper/single_dynamics_synth/pi2d_%d.png', k), 'png');    
    %     saveas(h, sprintf('/homes/49/edwin/matlab/combination/results/dynamicVB/paper/single_dynamics_synth/pi2d_%d.fig', k), 'fig');    
    end
    
    plot([0 0], [0 1], 'Color', 'black');    
    title('Evolution of E[\pi_j^k] for Deteriorating Agents');    
    grid on
end