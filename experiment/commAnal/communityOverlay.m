%Overlays two types of community: draws a bar graph showing the proportion
%of members of communities in g that are also members of each community in
%gTasks.

%For each pi-community, record the task community number of each member
gOverlay = cell(1, length(g));

figure;

Y = zeros(numel(g), numel(gTasks1));
X = 1:size(Ptasks1, 2);

for c=1:length(g)
    %find the task community for each member of a pi community
    piComm = g{c};
    
    if numel(piComm) < 20
        continue;
    end
        
    [maxVal, maxIdx] = max(Ptasks1(piComm,:), [], 2);
    
    gOverlay{c} = maxIdx;
    
    [N, X] = hist(gOverlay{c}, X);
    %N = N ./ sum(N);
    Y(c, X) = N;
end

Y = Y(sum(Y,2)>0, :);

Y = Y ./ repmat(sum(Y, 1), size(Y,1), 1);
Y = Y(:, sum(Y,1)>0);

bar3(Y');
figure;
bar(Y');