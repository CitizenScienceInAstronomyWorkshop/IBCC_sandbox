function [ Y ] = piSampleCommunityOverlay( gSam, gSamPi, Ppi, taskAgentMap )
%UNTITLED Summary of this function gSamoes here
%   Detailed explanation gSamoes here

%Overlays two types of community: draws a bar gSamraph showingSam the proportion
%of members of communities in gSam that are also members of each community in
%gSamTasks.

%For each pi-community, record the task community number of each member
gSamOverlay = cell(1, length(gSam));

figure;

Y = zeros(numel(gSam), size(Ppi,2));
X = 1:size(Ppi,2);

for c=1:length(gSam)
    %find the task community for each member of a pi community
    samComm = gSam{c};
    
    agents = taskAgentMap(samComm, :)~=0;
    
    [I,J] = find(agents); %just need to know J, the column no. i.e. the agSament index
        
    [maxVal, maxIdx] = max(Ppi(J, :), [], 2);
    
    gSamOverlay{c} = maxIdx';
    
    [N, X] = hist(gSamOverlay{c}, X');
    %N = N ./ sum(N);
    Y(c, :) = N;
end

Y = Y ./ repmat(sum(Y, 1), size(Y,1), 1);
Y = Y(:, sum(Y,1)>0);

bar3(Y');
figure;
bar(Y');


end

