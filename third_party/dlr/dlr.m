function os = dlr(X,T,PlotOut,Adapt);
%os = dlr(X,T,PlotOut,Adapt);
% X inputs, T targets : 0 or 1 {classes] or -1 [no label known]
% OPTIONAL PlotOut, Adapt [1,0], Default is 0 1
%os : output structure
%y : pr(class)
%w : weights
%S : sigma
%s : x'Sx
%q : weight diffusion variance
%a : latent variable
%ar : adapt rate ie mean trace (S)
  
  if (nargin == 2)
    PlotOut = 0;
    Adapt = 1;
  end;
  
  if (nargin < 2)
    error('dlr needs X,T as inputs');
  end;
  
  [L,D] = size(X);
%   if  (L<D)
%     X = X';
%     [L,D] = size(X);
%   end;
  
  if (length(T) ~= L)
    error('X,T must be same length');
  end;

  % normalise the entire data set and augment with column of ones (for
  % bias weights)
  X = normalis(X,X);
  X = [X , ones(L,1)];
  dim = size(X,2); % dimensionality of weights

  % now pre-allocate the variables
  os.x = zeros(L,dim);
  os.z = zeros(L,1);
  os.y = zeros(L,1);
  os.pc = zeros(L,1);
  os.yup = zeros(L,1);
  os.aup = zeros(L,1);
  os.sup = zeros(L,1);
  os.s = zeros(L,1);
%    os.S = zeros(L,dim,dim);
  os.ar = zeros(L,1);
  os.w = zeros(L,dim);
  os.q = zeros(L,1);
  os.a = zeros(L,1);
  os.v = zeros(L,1);
    
  % now read in everything and can start processing
  
  % first need to set up the priors
  I = eye(dim);
  Ir = diag([ones(1,dim-1),0]); % no noise on bias term latent updates if
                                % we use this
  S = I;
  F = I;
  w = zeros(dim,1); % last w is the bias
  q = 0;
  v=0;
  dq = 0;
  
  for t=1:L % scroll through the data
    x = X(t,:)'; % current input
    w = F*w;
    Sp = F*S*F' + q*I; % current adjusted sigma
    [y,a,s,kappa] = lat2out(w,Sp,x); % get output from latent space
    if (T(t) == -1) % no label is known
      z = y; % set target equal to expected outcome
      dq = y*(1-y);
    else
      z = T(t);
      dq = 0;
    end;
    u = y*(1-y);
    %Snew = Sp -( (y*(1-y))/(1 + y*(1-y)*s) )*(Sp*x)*(Sp*x)';
    %wnew = w + Snew*x*(z-y); % Newton step
    %wnew = w + (Sp/(1+y*(1-y)*s))*x*(z-y); % EKF step
    Snew = Sp - ((u * kappa^2) / (1 + u * s *	kappa^2)) * (Sp * x) * ...
        (Sp * x)';
    wnew = w + (kappa * Sp / (1 + u * s *kappa^2)) * x * (z - y); % EKF step
    S = Snew;
    w = wnew;
    [yup,aup,sup,kappaup] = lat2out(w,Snew,x);
    uup = yup*(1-yup);
    v = v + (uup-u);
    if  (Adapt == 1)
      q = max(uup - u, 0) + dq;
    end;

    %now get priors on non (-1) labels and use as decision threshold
    f = find(T ~= -1);
    classOneThresh = sum(T(f) == 1)/length(T(f));
    
    os.x(t,:) = x;
    os.z(t) = z;
    os.y(t) = y;
    os.pc(t) = (y>classOneThresh);
    os.yup(t) = yup;
    os.aup(t) = aup;
    os.sup(t) = sup;
    os.s(t) = s;
%      os.S(t,:,:) = S;
    os.ar(t) = mean(trace(S));
    os.w(t,:) = w;
    os.q(t) = q;
    os.a(t) = a;
    os.v(t) = v;
       
    if (PlotOut == 1)
      tmax = max(1,t-50);
      decbdy(os.x(tmax:t,1:2),os.pc(tmax:t),os.w(t,1:3),os.z(tmax:t)); %just plot projections onto first two dimensions
      axis([-3 3 -3 3]);
      drawnow;
    else
      if (rem(t,100)==0)
	fprintf('.');
      end;
    end;
  end;
  
  if (PlotOut == 0)
    fprintf('\n');
  end;
   
  
function [y,a,s,kappa] = lat2out(w,S,x)
    a = w'*x; % latent activation
    s = x'*S*x;
    y = moderate(a,s);
    kappa = (1 + (pi*(s)/8)).^(-0.5);
   