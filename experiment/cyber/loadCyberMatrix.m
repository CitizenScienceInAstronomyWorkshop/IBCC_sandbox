function [ cyMatrix, labels ] = loadCyberMatrix( filename, zeroScoreMap )
[fid, message] = fopen(filename);
cyRaw = textscan(fid, '%d %d %d %d %d %d %d %d %*s %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d', 'delimiter', ' :');

nRawRows = length(cyRaw{1});

labels = zeros(1, nRawRows);

maxAgentId = 0;
for c=9:length(cyRaw)
    maxCol = max(cyRaw{c});
    if maxCol > maxAgentId
        maxAgentId = maxCol;
    end
end
if nargin<2 || zeroScoreMap==0
    cyMatrix = sparse(double(maxAgentId), double(nRawRows));
else
    cyMatrix = zeros(double(maxAgentId), double(nRawRows)) + zeroScoreMap;
end
for r=1:nRawRows
    
    for j=1:8
        if cyRaw{j}(r) == 1
            labels(r) = j;
            break;
        end
    end
    for c=9:length(cyRaw)
        if mod(c-9, 2)==1 || isnan(cyRaw{c}(r)) || cyRaw{c}(r)==0
            continue;
        end
        cyMatrix(cyRaw{c}(r), r) = 1;
    end
end

fclose(fid);
end

