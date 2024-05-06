function ssr = loss_function_transformed(v2,X1,X0,Z1,Z0,w0)
  if nargin < 6
      w0 = [];
  end
  v = [1;exp(v2)]; % retransform v2
  D = diag(v);
  H = X0'*D*X0;
  H = (H + H')/2; % make sure that H is symmetric
  f = - X1'*D*X0;
  l = size(Z0,2);
  options = optimoptions('quadprog', 'Display', 'off');
  w = quadprog(H,f,[],[],ones(1,l),1,zeros(l,1),ones(l,1),w0,options);
  w = abs(w);
  e = Z1 - Z0*w;
  ssr = sum(e.^2);
