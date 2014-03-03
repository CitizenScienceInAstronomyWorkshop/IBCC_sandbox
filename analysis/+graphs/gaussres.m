function [px] = gaussres (x,mu,sigma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   [px] = gaussres (x,mu,sigma)
%
%   computes m-D probability density from x values given
%   mean mu and standard deviation sigma
%              1                                     T       -1
%   p(x)= ------------------------- exp (-0.5  (x-mu)   Sigma    (x-mu)  )
%         (2*pi)^(d/2) |Sigma|^0.5
%
%  e.g: [X,Y] = meshdom(-2:.2:2, -2:.2:2);
%       mu=[0;0];sigma=[1 0;0 1];
%       for i=1:21, for j=1:21, 
%             p(i,j)=gaussmd([ X(i,j) Y(i,j)],mu,sigma); 
%       end; end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d=size(sigma,1);
if (size(x,1)~=d)  x=x'; end;
if (size(mu,1)~=d)  mu=mu'; end;

ndim=size(x,1);
N=size(x,2);

z=(x-mu(:,ones(1,N)))';
IS=inv(sigma);
DS=det(sigma);

for l = 1:N,  
    px(l) = exp(-0.5*z(l,:)*IS*z(l,:)')  ./ sqrt((2*pi)^ndim*DS);
end;