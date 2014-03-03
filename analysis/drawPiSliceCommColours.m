function drawPiSliceCommColours(Pi, Kt, n, P, g, avgPis, minCommSize, chosenComms)
%drawPiEvolution Draw a 3-D graph showing the evolution of Pi over time.
%Works only when there are 3 possible scores (e.g. GZSN). A separate
%graph is drawn for each target label value.

%n should be the data point (in the score set format) for which we wish to
%draw a time-slice. For each base classifier k we use the last data point at
%which k made a classification, so k does not appear until 
%n=first_classification_k and does not move until n>next_classification_k
Alpha0 = [0.5 0.3 0.05; 0.18 0.36 0.41];
Pi0 = Alpha0 ./ repmat(sum(Alpha0,2), 1, size(Alpha0,2));

if nargin < 6
    minCommSize = 1;
end

minN = 0;

if size(Pi, 2) ~=3
    display('Pi Evolution is only implemented for confusion matrices with 3 possible scores.');
end

NKfiltered = [];
Kt = Kt(1:n);
Kunique = unique(Kt); %the ids of the base classifiers

if length(Kunique) == size(Pi, 3)
    NKfiltered = 1:size(Pi,3);
else
    for k=Kunique'
        Nk = find(Kt==k); %the data points corresponding to this base classifier
        NKfiltered = [NKfiltered Nk(end)];    
    end
end

h=figure('Position', [1 1 1000 700]); 

legendStrings = cell(size(P,2), 1);

ploth = zeros(1, size(P,2));

j = 1;
    
subplot(1, 1, 1);hold all;

X = [];%zeros(1, length(NKfiltered));
Y = [];%zeros(1, length(NKfiltered));
S = [];%zeros(1, length(NKfiltered)); %size of nodes
C = []; %colours based on community membership
Comms = [];
if nargin > 5 
    commCentres = zeros(2, length(g));
    for c=1:size(P, 2)

        membVals = P(:,c);
        [maxMembVal centreID] = max(membVals);

        [centreMaxMemb centreMaxComm] = max(P(centreID,:));
        if centreMaxComm ~= c || (nargin>7 && ~ismember(c, chosenComms))
            display(['skipping cluster ' num2str(c)]);
            continue;
        else
            display(['cluster ' num2str(c) ' has centre ID ' num2str(centreID)]);
        end

        Pij = avgPis(:,:,c);
        y = Pij(1);
        x = Pij(2)./sin(deg2rad(60)) + y./sin(deg2rad(60)).*sin(deg2rad(30));
        commCentres(:, c) = [x;y];
    end    
end

[membVals discreteComms] = max(P, [], 2);
for nk=NKfiltered

    Pij = Pi(j,:,nk);

    %skip those still at the priors as they are uninformative
    if sum(sum(Pij==Pi0(j,:,:)))>0
        display([num2str(nk) 'stuck at prior - skipping']);
        continue;
    end

    set(gca, 'FontSize', 18);

%         y = Pij(2)*sin(deg2rad(60));
%         x = (1-Pij(1)-Pij(2)) + y*cot(deg2rad(60));
    y = Pij(1);
    x = Pij(2)./sin(deg2rad(60)) + y./sin(deg2rad(60)).*sin(deg2rad(30));
    %         x = (1-Pij(2))./sin(deg2rad(60)) + y./sin(deg2rad(60)).*sin(deg2rad(30));

    x = round(x * 100)/100;
    y = round(y * 100)/100;

    [commVal comm] = max(P(nk,:));%Kt(nk)

    commAgents = find(discreteComms==comm);
    if numel(commAgents) < minCommSize
        continue;
    end

    if nargin > 7 && ~ismember(comm,chosenComms)
        continue;
    end
%         c = colours(comm, :);
    clashIdx = intersect(find(X==x), find(Y==y));
    colourClash = 0;
    for clIdx=clashIdx
        if Comms(clIdx)==comm
            colourClash = clIdx;
        end
    end
    if isempty(clashIdx) || colourClash==0
        X = [X x];
        Y = [Y y];
        S = [S 8];
%             C = [C c'];
        Comms = [Comms comm];
    else 
        S(colourClash) = (S(colourClash)^2 + 16)^0.5;
    end        
end
    
[S, idxs] = sort(S, 2, 'descend');
X = X(idxs);
Y = Y(idxs);

colours = colormap(hsv( length(unique(Comms)) ));   
for i=1:length(S)
%         ploth(j) = plot([X(i) commCentres(1, Comms(i))], [Y(i) commCentres(2, Comms(i))],...
      col= colours( find(unique(Comms)==Comms(i)), : );
      ploth(Comms(i)) = plot(X(i), Y(i),...
        '-o', 'Color', col', 'MarkerFaceColor', col', 'MarkerSize', S(i), 'MarkerEdgeColor', [0 0 0]);  
      legendStrings{Comms(i)} = ['Community ' num2str(Comms(i))];
%           text(X(i), Y(i), num2str(Comms(i)));
end

height = 1./sin(deg2rad(60));
plot([0 height height/2 0], [0, 0, 1 0], ':', 'color', 'black');   
textYOffset = -0.05;
textXOffset = -0.07;
h1 = text(textXOffset, textYOffset, 'score=3');
h2 = text(1+textXOffset, textYOffset, 'score=1');
h3 = text(height/2+textXOffset, 1.1+textYOffset, 'score=-1');        

ploth = ploth(find(ploth));
legendStringsCompact = cell(1,length(ploth));
i2=1;
for i=1:length(legendStrings)
    if ~isempty(legendStrings{i})
        legendStringsCompact{i2} = legendStrings{i};
        i2 = i2+1;
    end
end

legend(ploth, legendStringsCompact);        
set(gca, 'YTick', []);
set(gca, 'XTickLabel', []);

set(h1, 'FontSize', 18);
set(h2, 'FontSize', 18);
set(h3, 'FontSize', 18);
% xlabel('t');
title(['\pi^k_' num2str(j) ' vectors after ' num2str(n) ' data points']);
axis([0 height 0 height]);    

saveas(h, sprintf('/homes/49/edwin/matlab/combination/results/dynamicVB/paper/gzsn_singlefold/taskComm_%d.png', n), 'png');    
saveas(h, sprintf('/homes/49/edwin/matlab/combination/results/dynamicVB/paper/gzsn_singlefold/taskComm_%d.fig', n), 'fig');

end