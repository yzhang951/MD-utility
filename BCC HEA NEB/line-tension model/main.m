%在同一位置计算整条曲线

clear all
constants;
tau_xy =  0:4:2000; % MPa
num_stresses = length(tau_xy);
calc_times = 10000;
enthalpy = zeros(calc_times,length(tau_xy));
alpha0 = zeros(calc_times);
beta0 = zeros(calc_times);
%H_initial_arr = zeros(calc_times,length(tau_xy));

position = 0;
Xn = -N:1:N;
std = 0.0327;
para0 = [0.34, 4];     % initial value for [alpha, beta]

for j = 1:calc_times
    %tic
    z0_init = 0;
    pb = std*randn(2*N+1,1)+DeltaU;
    para0 = [0.34, 4];
    for i = 1:num_stresses
        [para,z0,H,H_diff,H_initial] = test_position(para0,z0_init,tau_xy(i),pb,0);
        if abs(imag(H))>1e-6 || real(H)<0
            break
        end
        if i == 1
            alpha0(j) = para(1);
            beta0(j) = para(2);
        end

        enthalpy(j,i) = real(H);
        %H_initial_arr(j,i) = real(H_initial);
        para0 = para; % use the saddle point of this step as initial value for the next calc
        z0_init = z0;
    end
    %toc
    if mod(j,500) == 0
        disp(['---- calc No. ', num2str(j),' finished'])
    end
    if mod(j,2000) == 0
        period_e = enthalpy(j-2000+1:j,:);
        filename = [num2str(j),'-enthalpy.mat'];
        save(filename,"period_e");
        disp(['---- calc No. ', num2str(j),' saved']);
    end
end

filename = 'homogenized-1022.mat';
save(filename);