function drawPiEvolution2scores(Alpha, Kt)
%Kt is the list of classifier IDs for each sample in the sample order, i.e. C{1} from the scoreset passed to the combiners
%drawPiEvolution Draw a 3-D graph showing the evolution of Pi over time.
%Works only when there are 3 possible scores (e.g. GZSN). A separate
%graph is drawn for each target label value.

minN = 0;

Kunique = unique(Kt); %the ids of the base classifiers


normTerm = repmat(sum(Alpha,2), [1 size(Alpha,2) 1]);
Pi = Alpha ./ normTerm;


if size(Pi, 2) ~=2
    display('Pi Evolution 2scores is only implemented for confusion matrices with 2 possible scores.');
end

Kfiltered = [];

for k=Kunique'
    
    Nk = find(Kt==k); %the data points corresponding to this base classifier
    
    nNk = length(Nk);
    if nNk < minN
        continue;
    end
    
    Kfiltered = [Kfiltered k];
end

display(num2str(length(Kfiltered)));

for k=Kfiltered
        Nk = find(Kt==k); %the data points corresponding to this base classifier

        display(['number of data points for ' num2str(k) ': ' num2str(length(Nk))]);
        h=figure('Position', [1 1 1000 700]); hold all

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
                Z(j, z) = z;
            end    
            w = 50;
            plot(Z(j,:), X(j,:), '-', 'linewidth', 2, 'color', colours{j});  
        end

    %     for j=1:size(Pi,1)
    %         idxs = (1:floor(length(X)/w)) .* w;
    %         idxs = [idxs length(X)];
    %         
    %         arrowCol = [0 0 0];        
    %         for a=idxs
    %             d=1;
    %                   
    %             if a-d < 1
    %                 continue;
    %             end
    %             
    %             while (X(j, a)-X(j, a-d))^2+(Y(j, a)-Y(j, a-d))^2 < 0.02^2
    %                 d=d+1;
    %                 if d==w || d > 5 || a-d < 1
    %                     break
    %                 end
    %             end
    %             
    %             if X(j, a)==X(j, a-d) && Y(j, a)==Y(j, a-d)
    %                 continue;
    %             end
    %             arrowCol(2) = 1-a/length(X);
    %             arrowCol(3) = 1-a/length(X);
    %             arrowh([X(j, a-d) X(j, a)], [Y(j, a-d) Y(j, a)], arrowCol, 150);
    %         end
    %     end

        textYOffset = -0.05;
        textXOffset = -0.07;
        h1 = text(0.5+textXOffset, textYOffset, 'score=0');
        h2 = text(0.5+textXOffset, 1-textYOffset, 'score=1');

        for j=1:size(Pi,1)
            legendStrings{j} = ['\pi^k_' num2str(j-1)];    
        end

        legend(legendStrings);        
        set(gca,'ZColor','w');
        set(gca,'ZTick', []);        
        set(gca, 'YTick', []);
    %     set(gca, 'XTickLabel', []);
        set(h1, 'FontSize', 18);
        set(h2, 'FontSize', 18);
        xlabel('sample no. i');
        title(['Evolution of E[\pi_j^k] Vectors for Base Classifier k=' num2str(k)]);
        axis([0 length(Nk) 0 1]);

    %     saveas(h, sprintf('/homes/49/edwin/matlab/combination/results/dynamicVB/paper/single_dynamics_synth/pi2d_%d.png', k), 'png');    
    %     saveas(h, sprintf('/homes/49/edwin/matlab/combination/results/dynamicVB/paper/single_dynamics_synth/pi2d_%d.fig', k), 'fig');    
end

end