function drawPiEvolution3d(Alpha, Kt, gSlices, filterMapSlices, slices)
%drawPiEvolution Draw a 3-D graph showing the evolution of Pi over time.
%Works only when there are 3 possible scores (e.g. GZSN). A separate
%graph is drawn for each target label value.

normTerm = sum(Alpha, 2);
normTerm = repmat(normTerm, 1, size(Alpha,2));
Pi = Alpha ./ normTerm;

minN = 1000;

Kunique = unique(Kt); %the ids of the base classifiers

if size(Pi, 2) ~=3
    display('Pi Evolution is only implemented for confusion matrices with 3 possible scores.');
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

    h=figure('Position', [1 1 1600 500]);
    subplot(1, 2, 1);
    w = 50;    
    for z=1:length(Nk)

        if mod(z, w) ~=0
            continue;
        end
        height = 1./sin(deg2rad(60));
        
        plot3([z z z z], [0 height/2 height 0], [0 height 0 0], '-', 'color', [0.9 0.9 0.9]); hold on
    end 
    
    legendStrings = cell(size(Pi,1), 1);     
   
    colours = {[0,0.3,0.3], [1,0,0]};
    
    X = zeros(size(Pi,1), length(Nk));
    Y = zeros(size(Pi,1), length(Nk));
    Z = zeros(size(Pi,1), length(Nk));    
    
    hline = zeros(1,size(Pi,1));
       
    for j=1:size(Pi,1)

        Pij = Pi(j,:,:);

        set(gca, 'FontSize', 14);

        for n=Nk'
            Pijt = Pij(1, :, n);
            y = Pijt(1);
            x = Pijt(3)./sin(deg2rad(60)) + y./sin(deg2rad(60)).*sin(deg2rad(30));
            z = find(Nk==n);

            X(j, z) = x;
            Y(j, z) = y;
            Z(j, z) = z;
        end    
%         hline(j)=plot3(Z(j,:), X(j,:), Y(j,:), '-', 'linewidth', 2, 'color', colours{j}); hold all  
    end
          
    textOffset = -170;    
    h1 = text(textOffset, 0,  0, 'score=1');
    h2 = text(textOffset, 1.1, 0, 'score=3');
    h3 = text(textOffset, 0.5, 1.1, 'score=-1');
%     for j=1:size(Pi,1)
%         legendStrings{j} = ['\pi^k_' num2str(j-1)];    
%     end
      
    hold on;  
    
    plot3([0 0 0 0 ], [0 height/2 height 0], [0 height 0 0], '-', 'color', [0 0 0]); hold on
    plot3([length(Nk) length(Nk) length(Nk) length(Nk)], [0 height/2 height 0], [0 height 0 0], '-', 'color', [0 0 0]); hold on
    plot3([0 length(Nk)], [height/2 height/2], [height height], '-', 'color', [0 0 0]); hold on
    plot3([0 length(Nk)], [height height], [0 0], '-', 'color', [0 0 0]); hold on
    plot3([0 length(Nk)], [0 0], [0 0], '-', 'color', [0 0 0]); hold on

    set(gca,'ZTick', []);        
    set(gca, 'YTick', []);
    set(gca, 'XTickLabel', []);
    set(h1, 'FontSize', 14);
    set(h2, 'FontSize', 14);
    set(h3, 'FontSize', 14);
    xlabel('timestep \tau');
    title(['Changes in Confusion Matrix of Base Classifier ' num2str(k)]);
    axis([0 length(Nk) 0 height 0 1]);
    
    legendStrings{1} = 'not supernova';
    legendStrings{2} = 'supernova';
    
    for s=1:numel(slices)
        Nk = find(Kt==k);
        if k>length(filterMapSlices{s})
            sliceId = 0
        else
            sliceId = filterMapSlices{s}(k);
        end
        if sliceId > 0
            x = max(find(Nk<slices(s)));
            
            Pit = Pi(1,:,Nk(x));
            
            z = Pit(1);
            y = Pit(3)./sin(deg2rad(60)) + z./sin(deg2rad(60)).*sin(deg2rad(30));            

            Pit = Pi(2,:,Nk(x));
            
            z2 = Pit(1);
            y2 = Pit(3)./sin(deg2rad(60)) + z2./sin(deg2rad(60)).*sin(deg2rad(30));            
            plot3([x x],[y y2],[z z2],'-*', 'MarkerSize', 12);
            
            comm = 0;
            for c=1:numel(gSlices{s})
                if ~isempty(find(gSlices{s}{c}==sliceId))
                    comm = c;
                    break;
                end
            end
            if comm>0
                text(x,y,z-0.05,num2str(slices(s)), 'FontSize', 16);
            end
        end
    end
        
    for j=1:size(Pi,1)
        hline(j)=plot3(Z(j,1), X(j,1), Y(j,1), '*', 'linewidth', 2, 'color', colours{j}); hold all  
    end
    
    hLeg = legend([hline(1) hline(2)], legendStrings);    
    set(hLeg, 'Location', 'northeast');
    
    clear M;
    
    mkdir('/homes/49/edwin/results/ahm12_demo/toblerone_playground/', num2str(k));
        
    display(num2str(Nk));
    for p=1:length(Nk)
        subplot(1, 2, 1);
        if mod(p,10)~=0
            continue;
        end
        
        for j=1:size(Pi,1)
            delete(hline(j));
            hline(j)=plot3(Z(j,1:p), X(j,1:p), Y(j,1:p), '-', 'linewidth', 3, 'color', colours{j}); hold all  
        end
        
        %the 2d version
        subplot(1, 2, 2);
        hold off

        plot([0 0.5*height height 0], [0 height 0 0], '-', 'linewidth', 1, 'color', 'black'); hold all
        for j=1:size(Pi,1)
            plot(height-X(j,p), Y(j,p), '.', 'linewidth', 3, 'MarkerSize', 50, 'color', colours{j}); hold all  
        end       
        title('Cross-sectional snapshot view of current time-step');    
        textOffset = 0.2;    
        h1 = text(height+textOffset/4, textOffset/4, 'score=1');
        h2 = text(0.5, height+textOffset/4, 'score=-1');
        h3 = text(-textOffset, textOffset/4, 'score=3');
        set(h1, 'FontSize', 14);
        set(h2, 'FontSize', 14);
        set(h3, 'FontSize', 14);         
        
        M(p) = getframe(h);
        imwrite(M(p).cdata, sprintf('/homes/49/edwin/results/ahm12_demo/toblerone_playground/%d/%.4d.png', k, p))
%         saveas(h, sprintf('/homes/49/edwin/results/ahm12_demo/toblerone_playground/%d/%.4d.png', k, p), 'png');   
        display(num2str(p)); %frames completed
    end 
    
%     movie2avi(M, sprintf('/homes/49/edwin/results/ahm12_demo/toblerone_playground/%d/3dMovie_%d',k));
%     saveas(h, sprintf('/homes/49/edwin/results/ahm12_demo/toblerone_playground/pi_%d.png', k), 'png');    
%     saveas(h, sprintf('/homes/49/edwin/results/ahm12_demo/toblerone_playground/pi_%d.fig', k), 'fig'); 
%     saveas(h, sprintf('/homes/49/edwin/results/ahm12_demo/toblerone_playground/pi_%d.eps', k), 'epsc'); 
end

end