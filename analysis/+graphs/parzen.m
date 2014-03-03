%[m,C,p] = parzen(x);

function [m,C,p] = parzen(x)

TAU = 4;

[K,D] = size(x);
m = x;
p = ones(K,1)/K;

for k = 1:K
  temp = (m-ones(K,1)*m(k,:));
  temp = temp.^2;
  temp = sum(temp);
  %temp = sort(sum(temp'));
  %temp=temp(2);
  C(k,:)=ones(1,D).*temp;
end;

C = C*TAU;