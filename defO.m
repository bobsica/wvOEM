function O = defOreal

[O] = oem;
O.A = true;
O.cost = true;
O.dx = true;
O.e = true;
O.eo = true;
O.es = true;
O.ex = true;
O.G = true;
O.J = true;
O.jexact = true;
O.jfast = false;
O.linear = false;
%
if ~O.linear 
  O.itermethod = 'ML'; % ?GN? Gauss-Newton ?ML? or ?LM? for Marquardt-Levenberg
  O.stop_dx = .05; % .1 .005;.1; %1; %2; %.5, 1,50;
  O.maxiter = 10; % 8, 12
  O.ga_factor_not_ok = 10;
  O.ga_factor_ok = 10;
  O.ga_max = realmax; %1e15,1e12,1e10;%1e8,1e6
  O.ga_start = 100; % 100
end
O.S = true;
O.So = true;
O.Ss = true;
O.sxnorm = true;
O.yf = true;
O.Xiter = true; 

return
