function [ data ] = selectPortionUnknown( nKnown1, nRatioUnknown, data, knownAssets)
%SELECTPORTIONUNKNOWN Selects a portion of the unknown labels

%     nKnown2Desired = nKnown1*nRatioKnown2;
    nUnknown = nKnown1*nRatioUnknown; %have more unknowns than known2s
%     agents = zeros(nUnknown+nKnown1+nKnown2Desired, 1);
%     assets = zeros(nUnknown+nKnown1+nKnown2Desired, 1);
%     scores = zeros(nUnknown+nKnown1+nKnown2Desired, 1);  
%     
%     u=0; k=0; k2=0;
    
    unknownAssets = find(knownAssets==0);
    unknownAssets = unknownAssets(1:nUnknown);
    knownAssets = find(knownAssets~=0);
    assetsToKeep = [unknownAssets knownAssets];
    idxsToKeep = ismember(data{2}, assetsToKeep);
    agents = data{1}(idxsToKeep);
    assets = data{2}(idxsToKeep);
    scores = data{3}(idxsToKeep);
%     
%     for i=1:length(data{2})        
%         if sum(knownAssets1==data{2}(i))>0
%             agents(u+k+k2+1) = data{1}(i);
%             assets(u+k+k2+1) = data{2}(i);
%             scores(u+k+k2+1) = data{3}(i);
%             k = k+1;
%         elseif sum(knownAssets2==data{2}(i))>0
%             if k2<k*nRatioKnown2
%                 agents(u+k+k2+1) = data{1}(i);
%                 assets(u+k+k2+1) = data{2}(i);
%                 scores(u+k+k2+1) = data{3}(i);
%                 k2 = k2+1;
%             else
%                 continue
%             end
%         elseif u<k*nRatioUnknown
%             agents(u+k+k2+1) = data{1}(i);
%             assets(u+k+k2+1) = data{2}(i);
%             scores(u+k+k2+1) = data{3}(i);
%             u = u+1;
%         end        
%     end
%     
%     display(num2str(u));
    
    data = {double(agents) double(assets) double(scores)};
end

