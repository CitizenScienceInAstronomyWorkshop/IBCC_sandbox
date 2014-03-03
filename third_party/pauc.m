function [p] = pauc (auc1,auc2,n0,n1)

% function [p] = pauc (auc1,auc2,n0,n1)
% Return probability of auc2 > auc1 if each is an AUC statistic
% See Hanley and McNeil

% auc1 is always set to 0.5, change func def later

%auc1=0.5;

% Do this, as we dont know a priori which variable will be greater
auc2=max(auc2,1-auc2);

q1=auc2/(2-auc2);
q2=2*auc2^2/(1+auc2);
se=(auc2*(1-auc2)+(n1-1)*(q1-auc2^2)+(n0-1)*(q2-auc2^2))/(n0*n1);
se=sqrt(se);

z=(auc2-auc1)/se;
p=pnorm(z);
