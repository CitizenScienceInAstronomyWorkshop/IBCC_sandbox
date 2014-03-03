function graphCombinationInProgress( testIdxs, results, testLabels, foldData)
%Displays responses from base classifiers alongside combined responses
    count = 0;
    countPos = 0;
    countNeg = 0;
    
    for i=randi(length(testIdxs), 1, length(testIdxs))
        
        if abs(results(2,i)-(testLabels(i)-1)) > 0.5
            continue;
        end
        
        if testLabels(i)==2 && countPos-countNeg>5
            continue;
        end
        
        if testLabels(i)==1 && countNeg-countPos>5
            continue;
        end
        
        if testLabels(i)==2
            countPos = countPos+1;
        end
        
        if testLabels(i)==1
            countNeg = countNeg+1;
        end
        
        count = count+1;
        
        if count > 30
            break
        end
        
        h = figure('Position', [1 1 1600 500]);                
        
        subplot(1,2,1);
                
        foldDataIdxs = find(foldData{2}==testIdxs(i));
        nResp = length(foldDataIdxs);
%         plot([0 1], [foldData{3}(foldDataIdxs) foldData{3}(foldDataIdxs)], '.-', 'LineWidth', 10, 'MarkerSize', 5);

        hold all

        for r=1:nResp
            bar(r, foldData{3}(foldDataIdxs(r)));
        end
        axis([0 nResp+1 -1 3]);
        set(gca, 'FontSize', 14);
        
        title('Volunteer Responses for a sample candidate from GZ Supernovae'); 
        xlabel('Index of Responders');
        ylabel('Scores');
        
        subplot(1,2,2);
        hold all
        bar(1, (results(1,i)+1)/4, 'y');
        bar(2, results(2,i), 'b');
        bar(3, testLabels(i)-1, 'r');
              
        axis([0 4 0 1]);
        set(gca, 'FontSize', 14);
        set(gca,'XTickLabel',{'Mean', 'IBCC', 'True Label'});
        set(gca,'XTick', [1 2 3]);
        ylabel('Probability of supernova');
        title('Combined Decisions and True Label');
        
        f = getframe(h);
        imwrite(f.cdata, sprintf('/homes/49/edwin/results/ahm12_demo/scoresScrollThrough/%.4d.png', count))        
        close all
    end

end

