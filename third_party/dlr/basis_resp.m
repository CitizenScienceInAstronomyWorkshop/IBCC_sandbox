function [B,m] = basis_resp(x, Nh, m, g_or_tps, verbose, add_linear, g_var);
% [B,m] = basis_resp(x, Nh, OPT m, OPT g_or_tps, OPT verbose, OPT add_linear, OPT g_var);
% makes a basis and then returns the Hilbert expansion
% B = [x, phi(x)] for a set of Nh basis responses
% a vector of ones is then appended in mapping routines to make
% Bfinal = [x, phi(x), 1] as universal Hilbert projection
% (c) Stephen Roberts (2007-)
  
  if (nargin < 4)
      g_or_tps = 'g';
  end;
  if (nargin < 5)
      verbose = 1;
  end;
  if (nargin < 6)
      add_linear = 1;
  end;
     
  if (g_or_tps ~= 'g') & (g_or_tps ~= 'tps')
      error('wrong basis function string, g or tps');
  end;
  
  if (verbose == 1)
      if (g_or_tps == 'g')
          disp('gaussian');
      else 
          disp('thin plate spline');
      end;
  end;
    
  [nts,dts] = size(x);
  B = [];
  
  if (Nh > 0)
      if (nargin < 3) | (isempty(m) == 1);
          for n=1:Nh
              c = ceil(rand*nts);
              m(n,:) = x(c,:) + 0.1*(randn(1,dts));
          end;
      else
          if (verbose==1) disp('using input m'); end;
      end;

      if (nargin < 7)
          sig2 = 100*max(diag(cov(m)));
      else
          sig2 = g_var;
      end;
  
      for n=1:Nh
          z =  ((x-ones(nts,1)*m(n,:))'.^2);
          if (dts > 1)
              z = sum(z);
          else
              z = z';
          end;
          if (g_or_tps == 'tps')
              phi(:,n) = (0.5*z.*log(z+eps))';
          else
              phi(:,n) = exp(-0.5*z/sig2);
          end;
      end;
      if (add_linear == 1)
          B = [x,phi];
      else
          B = phi;
      end;
  else % Nh = 0, so no basis expansion. Default add in linear
      B = x;
  end;