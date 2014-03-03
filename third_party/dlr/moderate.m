% [y] = moderate(a,v)
%
% a	activations
% v          variances
% y          resultant logistic outputs
% (c) Stephen Roberts (2005-)

function [y] = moderate(a,v)

[L,D] = size(a);

if (D == 1)
    a = [a, 0];
    v = [v, v];
end;

if (size(a) ~= size(v))
  error('a & v must be of the same size!');
  return;
end;

k = (1 + (pi*(v)/8)).^(-0.5);
y = softmax(k.*a);

if (D == 1)
    y = y(:,1);
end;

