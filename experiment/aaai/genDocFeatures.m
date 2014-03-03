
nF = 250;
Xmat = zeros(nDocs, nF);

nTopics = length(chosenIdx);
nFPerTopic = 50;

nSpikes = 3;
nSpikers = 15;

for d=1:nDocs
    
    topic = find(qRels(selectedDocs(d),:));
   
    topic = find(ismember(topic, chosenIdx));
    
    for t=1:nTopics
        start = (t-1)*nFPerTopic +1;
        last = t*nFPerTopic;
        if t==topic 
            
            spike = randi(nSpikers,nSpikes,1);
            spike2 = randi(nSpikers,nSpikes,1);
            spike3 = randi(nSpikers,nSpikes,1);
            spike4 = randi(nSpikers,nSpikes,1);
            
            
            Xmat(d, start:last) = betarnd(1,7,1,last-start+1);
            Xmat(d, start+spike) = betarnd(5,1,1,nSpikes);
            Xmat(d, start+nSpikers+spike2) = betarnd(3,2,1,nSpikes);
            Xmat(d, start+nSpikers+nSpikers+spike3) = betarnd(1,2,1,nSpikes);
            Xmat(d, start+nSpikers+nSpikers+nSpikers+spike4) = betarnd(3,2,1,nSpikes);

        else
            Xmat(d, start:last) = betarnd(1,7,1,last-start+1);
        end
    end
    
    if isempty(topic)
        spike = randi(nFPerTopic,nSpikes,1);
        
        Xmat(d, last+1:last+nFPerTopic) = betarnd(1,7,1,nFPerTopic);
        Xmat(d, last+1+spike) = betarnd(5,1,1,nSpikes);
                    
        spike2 = randi(nSpikers,nSpikes,1);
        Xmat(d, start+nSpikers+spike2) = betarnd(3,2,1,nSpikes);
        
        spike3 = randi(nSpikers,nSpikes,1);
        Xmat(d, start+nSpikers+nSpikers+spike3) = betarnd(1,2,1,nSpikes);
        
        spike4 = randi(nSpikers,nSpikes,1);
        Xmat(d, start+nSpikers+spike2+nSpikers+nSpikers) = betarnd(3,2,1,nSpikes);
    else
        Xmat(d, last+1:last+nFPerTopic) = betarnd(1,7,1,nFPerTopic);
    end
    
     Xmat(d, last+nFPerTopic+1:end) = betarnd(1,1,1,nF-last-nFPerTopic);
end

[firstCol secCol] = ind2sub(size(Xmat), 1:numel(Xmat));

Xsynth = zeros(size(firstCol,2), 3);
Xsynth(:,1) = selectedDocs(firstCol)'-1;
Xsynth(:,2) = secCol'-1;
Xsynth(:,3) = reshape(Xmat, numel(Xmat), 1);

Xsynth = [[max(selectedDocs) 0 0; nF 0 0]; Xsynth];