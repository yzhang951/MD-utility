function z0 = calc_initial(z0_init,tau_xy,pb)

constants;
%tau_xy = 100;
%std = 0.01 ;
%DeltaU = 0.03487;
%pb = std*randn(2*N+1,1)+DeltaU;
fun = @(z0)initial_state(z0,tau_xy,pb);
% algorithrms:
%   trust-region
%   levenberg-marquardt
%   trust-region-dogleg
%z0_init = 0;
options = optimoptions('fsolve','Algorithm','trust-region-dogleg','FunctionTolerance',1e-20,'OptimalityTolerance',1e-20,'StepTolerance',1e-20,'Display','off');
z0=fsolve(fun,z0_init, options);

%z = 1/2/pi*(asin(tau_xy/converter*b^2*Lp/DeltaU/pi))

end
%%