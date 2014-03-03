% [y] = softmax(a)
%
% a	activations
% y resultant logistic outputs
% (c) Stephen Roberts, (2000-)

function [y] = softmax(a)

[N,D] = size(a);
if (D==1)
    a = [a, 0];
end;
e = exp(a);
s = sum(e')';
y = e ./ (s*ones(1,D));

if (D==1)
    y = y(:,1);
end;
