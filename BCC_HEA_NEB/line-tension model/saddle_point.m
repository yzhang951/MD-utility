function [para,z0,H,H_diff,H_initial] = saddle_point(para0,z0_init,tau_xy,pb,position)
constants;

z0 = calc_initial(z0_init,tau_xy,pb);
H_initial = calc_H([0,0],tau_xy,pb,z0,0);
fun = @(para)diff_H(para,tau_xy,pb,z0,position);

options = optimoptions('fsolve','Algorithm','trust-region','FunctionTolerance',1e-20,'OptimalityTolerance',1e-20,'StepTolerance',1e-20,'Display','off');
para=fsolve(fun,para0, options);

H = calc_H(para,tau_xy,pb,z0,position)-H_initial;
H_diff = diff_H(para,tau_xy,pb,z0,position);

end