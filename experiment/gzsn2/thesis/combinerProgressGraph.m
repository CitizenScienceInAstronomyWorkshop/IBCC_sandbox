

logs = {};
SD = progressMonitor.calculateDeviation();
logsSD = {};
%SKIP any empty error logs for non-iterative methods
for c=1:length(progressMonitor.error)
    if ~isempty(progressMonitor.error{c})
       logs{length(logs)+1} = progressMonitor.error{c};
       logsSD{length(logsSD)+1} = SD{c};
    end
end

nFigs = size(logs{1},1);
for f=1:nFigs
    
    figure;
    
    for c=1:length(logs) 
        error = logs{c};
        e = error(f, :);
        lowerBar = e - logsSD{c}(f,:);
        upperBar = e + logsSD{c}(f,:); 
        lineHandle = plot(e, 'linewidth', 2); hold on;
        col = get(lineHandle, 'color');
        barIdxs = c*5:20:size(lowerBar,2);
        for i=barIdxs
            plot([i i], [lowerBar(i) upperBar(i)], 'x-', 'linewidth', 1, 'color', col); hold on
        end
        
        display('Also adding in a shortened version so you can do a close-up of the first few iterations');
%         if length(e)>100
%             hold all
%             plot(e(1:100),'linewidth',2);
%         end
        hold all;
    end

    fs = 24;
    set(gca, 'FontSize', fs);
    if length(logs)==2
        legend({'IBCC-VB', 'IBCC-Gibbs'});
    elseif length(logs)>2
        legend({'IBCC-VB-SubClass', 'IBCC-VB', 'IBCC-Gibbs'});
    end
    xlabel('No. Iterations');
    
    if f==nFigs
        ylabel('Entropy');
        title('Decreasing Uncertainty in Target Labels');
    else
        ylabel('AUC');%'Total absolute error over 1000 data points');
        title(['Class ' num2str(f)]);%Absolute error against number of iterations of IBCC methods');
    end
    
    set(gca, 'fontsize', fs) 
    
    grid on
end