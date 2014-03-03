nClasses = 8;
nRequests = 1000;

figure;

%limited to 300 because the later points are all basically the same


% data = lam_uncert_onceonly(1,1:nRequests)'; 
% xpoints = [1:nRequests].*3;
% plot(xpoints, data);    
% hold all
% 
% data = lam_uncert_5timesfixed(1,1:nRequests)'; 
% xpoints = [1:nRequests].*3;
% plot(xpoints, data);    
% hold all
% 
% data = lam_uncert_unrestricted(1,1:nRequests)'; 
% xpoints = [1:nRequests].*3;
% plot(xpoints, data);    
% hold all


% data = fpr_comb(1,1:nRequests)';
% xpoints = [1:nRequests].*3;
% plot(xpoints, data);
% hold all

data = lam_comb2(1,1:nRequests)';
xpoints = [1:nRequests].*3;
plot(xpoints, data);
hold all

data = lam_knockon(1,1:nRequests)';
xpoints = [1:nRequests].*3;
plot(xpoints, data);
hold all

data = lam_uncert(1,1:nRequests)'; 
xpoints = [1:nRequests].*3;
plot(xpoints, data);    
hold all    

% data = fpr_rand(1,1:nRequests)'; 
% xpoints = [1:nRequests].*3;
% plot(xpoints, data);    
% hold all
% 
% data = fpr_rand2(1,1:nRequests)'; 
% xpoints = [1:nRequests].*3;
% plot(xpoints, data);    
% hold all
% 
data = lam_rand3(1,1:nRequests)'; 
xpoints = [1:nRequests].*3;
plot(xpoints, data);    
hold all

% legendStrings = {'Approx 3-step lookahead v1', 'Approx 3-step lookahead v2', '1-step IG', 'Most uncertain',...
%     'Random, 1st run', 'Random, 2nd run', 'Random, 3rd run'};
legendStrings = {'Approx 3-step lookahead v2', '1-step IG', 'Most uncertain',...
    'Random, 1st run',};
% legendStrings = {'once only', 'five times', 'unrestricted'};
legend(legendStrings);

xlabel('number of labels requested');
ylabel('LAM');
title('Second Subset of Cyber Data');