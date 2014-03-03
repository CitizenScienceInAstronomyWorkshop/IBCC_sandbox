function h = drawPiSlice(Pi, g, P, j, gIdx)
%drawPiEvolution Draw a 3-D graph showing the evolution of Pi over time.
%Works only when there are 3 possible scores (e.g. GZSN). A separate
%graph is drawn for each target label value.

%n should be the data point (in the score set format) for which we wish to
%draw a time-slice. For each base classifier k we use the last data point at
%which k made a classification, so k does not appear until 
%n=first_classification_k and does not move until n>next_classification_k

minN = 0;

if size(Pi, 2) ~=3
    display('Pi Evolution is only implemented for confusion matrices with 3 possible scores.');
end

h=figure('Position', [1 1 1000 900]); 
hold all;
legendStrings = cell(length(g), 1);
% colours = {[0,0,0], [1,0,0]};
colours = colormap(hsv(size(Pi,3)));

ploth = zeros(1, size(Pi,1));

style = {' o', ' o'};


    
    X = [];%zeros(1, length(NKfiltered));
    Y = [];%zeros(1, length(NKfiltered));
    S = [];%zeros(1, length(NKfiltered)); %size of nodes
    C = [];
    Comms = [];

    for comm=1:size(Pi,3)
    
        Pij = Pi(j,:,comm);

        set(gca, 'FontSize', 22);

%         y = Pij(2)*sin(deg2rad(60));
%         x = (1-Pij(1)-Pij(2)) + y*cot(deg2rad(60));
        y = Pij(1);
%         x = (1-Pij(2))./sin(deg2rad(60)) + y./sin(deg2rad(60)).*sin(deg2rad(30));
        x = Pij(2)./sin(deg2rad(60)) + y./sin(deg2rad(60)).*sin(deg2rad(30));
        
        x = round(x * 100)/100;
        y = round(y * 100)/100;
        
        clashIdx = intersect(find(X==x), find(Y==y));
        c= colours(comm, :);
%         if isempty(clashIdx)
            X = [X x];
            Y = [Y y];
            S = [S 3.5*length(g{gIdx(comm)})^0.5]; %use three for when there are lots of overlapping large communities
            C = [C c];
            Comms = [Comms comm];
%         else 
%             X = [X x];
%             Y = [Y y];
%             S = [S S(clashIdx(end))+2];
%             C = [C c];            
%             S(clashIdx) = (S(clashIdx)^2 + 16)^0.5;
%         end        
    end
    
    [S, idxs] = sort(S, 2, 'descend');
    X = X(idxs);
    Y = Y(idxs);
    Comms = Comms(idxs);
    
    for i=1:length(S)
        ploth(j) = plot(X(i), Y(i), style{j}, 'MarkerFaceColor', colours(i, :), 'MarkerSize', S(i), 'MarkerEdgeColor', [0 0 0]);  
%         text(X(i), Y(i), num2str(Comms(i)));
    end

    
height = 1./sin(deg2rad(60));
plot([0 height height/2 0], [0, 0, 1 0], ':', 'color', 'black');   
textYOffset = -0.05;
textXOffset = -0.07;
h1 = text(textXOffset, textYOffset, 'score=3');
h2 = text(1.15+textXOffset, textYOffset, 'score=1');
h3 = text(height/2+textXOffset, 1.1+textYOffset, 'score=-1');    
% 
% for j=1:length(g)
% %     legendStrings{j} = ['\pi^k_' num2str(j)];    
%     legendStrings{j} = ['Community ' num2str(Comms(j))];
% end

% legend(legendStrings);        
set(gca, 'YTick', []);
set(gca, 'XTickLabel', []);

set(h1, 'FontSize', 22);
set(h2, 'FontSize', 22);
set(h3, 'FontSize', 22);
% xlabel('t');
if j==1
    title('Not Supernova');
elseif j==2
    title('Supernova');
else
    title(['Mean \pi_j^k Vectors']);% after ' num2str(n) ' data points']);
end
axis([0 height 0 height]);

% saveas(h, sprintf('/homes/49/edwin/matlab/combination/results/dynamicVB/paper/gzsn_singlefold/slice_%d.png', n), 'png');    
% saveas(h, sprintf('/homes/49/edwin/matlab/combination/results/dynamicVB/paper/gzsn_singlefold/slice_%d.fig', n), 'fig');

end