classdef IbccVbSeqPriorUpdates < combiners.bcc.IbccVb
    %IBCCVBSEQPRIORUPDATES Modifies combineScoreSet so that sequential
    %updates us a prior generated from the previous updates, rather than
    %updating the whole batch.
    
    properties
    end
    
    methods
        function [ET, nIt, EPi, EAlpha] = combineScoreSet(obj, C, Tvec, Alpha)
            nSamples = length(C{2});
            nAssets = length(Tvec);
            nIt = 0;
            display(['IBCCVB: number of samples to process: ' num2str(nSamples)]);
            
            if obj.sequential
                firstIt = obj.nBootstrap;
                if firstIt > nSamples
                    firstIt = nSamples;
                end
                if obj.seqPriorUpdates
                    start = 1;               
                end
            else
                firstIt = nSamples;                            
            end

            for i=firstIt:nSamples
                display(['sample no: ' num2str(i) ', no. VB iterations: ' num2str(nIt)]);
                
                if obj.seqPriorUpdates
                    [nIt, start] = obj.updateScoreSetSPU(start, i, C, post, Tvec);
                    continue;
                end            
            end
            
            EPi = exp(lnPi);
        end 
        
        %sequential IBCC with prior updates rather than continuing the iterations
        function [nIt, start] = updateScoreSetSPU(obj, start, i, C, post, T)

            %won't work with scoreset
            C_i = C(:, start:i);
            post_i = post(:, start:i);
            T_i = T(start:i);

            %change T into the matrix format for multi-class problems
            %update the initial values
            [Alpha, lnPi, lnP, ET] = obj.initVariables(T_i, C_i, post_i, Alpha);

            if obj.useKnownsOnly
                C_i = C;
                labelsKnown = find(T_i(C_i{2})~=0);
                C_i{1} = C_i{1}(labelsKnown);
                C_i{2} = C_i{2}(labelsKnown);
                C_i{3} = C_i{3}(labelsKnown);
            end                

            converged = false;
            threshold = obj.convThreshold;

            nItOld = nIt;

            while ~converged                
                oldET = ET;
                nIt = nIt + 1;
                [Alpha, lnPi, lnP, ET, pT] = obj.iterate(Alpha, lnPi, lnP, T_i, C_i);

                if sum(abs(ET(2,:)-oldET(2,:))) < obj.convThreshold || ...
                        (nIt > nItOld+5 && i <= 0.5 * length(T))
                    ET = oldET;
                    cIt = cIt + 1;
                else
                    cIt = 0;
                end
                if cIt > obj.convIt  || nIt > obj.maxIt
                    converged = true;
                end
            end

            Counts = obj.voteCounts(C_i, ET);                    
            for j=1:obj.nClasses
                Alpha(j,:,:) = (Alpha(j,:,:) + Counts{j}); % .* obj.discountAlpha;
            end

            obj.nu = obj.nu + sum(ET, 2)';

            %need to change this for >2 classes as it only gives prob. of T==2
            obj.combinedPost(1, start:max(C_i{2})) = ET(2, :);

            %next time do not include the first known targets
            start = i+1;                 
        end     
        

        function combinedPost = combineCluster(obj, post, Tvec)          
            obj.nKnown = round(obj.nKnown);
            C = obj.prepareC(post);
            Alpha = obj.initVariables(Tvec, C);
            
            [ET, nIt, EPi, EAlpha] = obj.combineScoreSet(C, Tvec, Alpha);              
            combinedPost = obj.combinedPost;

            nItFile = sprintf('%s%s', obj.confSave, 'nIts_vb');
            if exist(nItFile, 'file');
                prevNIts = dlmread(nItFile);
                prevNIts = [prevNIts nIt];
            else
                prevNIts = nIt;
            end
%             dlmwrite(nItFile, prevNIts);
            
            display([num2str(nIt) ' iterations for IBCCVB. sequential=' num2str(obj.sequential) ' and prior updates=' num2str(obj.seqPriorUpdates)]);
            
            obj.Alpha = obj.Alpha + EAlpha;
        end
         
    end
    
end

