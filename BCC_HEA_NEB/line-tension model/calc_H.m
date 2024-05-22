function H=calc_H(para,tau_xy,pb,z0,position)
    % para(1) -> alpha
    % para(2) -> beta
    constants;
    tau_xy = tau_xy / converter; % convert MPa to ev/A^3
    alpha = para(1);
    beta = para(2);
    H = 0;
    X = (-N):1:N; % position along dislocation line
    %z0 = 0;%1/2/pi*(asin(tau_xy*b^2*Lp/DeltaU/pi));
    Z = zeros(1,2*N+1);
    for n = 1:(2*N+1)
        Z(n) = (0.5 * tanh(alpha * (X(n)-position + beta)) - 0.5 * tanh(alpha * (X(n)-position - beta)))*(1-z0)+z0;
    end
    
    for n = 1:(2*N+1)
        if n ~= 2*N+1                % periodic BC
            subs = Z(n) - Z(n+1);
        else
            subs = Z(n) - Z(1);
        end

        Up = pb(n)/2 * (1 - cos(2*pi*Z(n)/period));
        load = tau_xy*b^2*Lp * Z(n);
        lt = Gamma*Lp^2/2/b * (subs)^2;

        H = H + Up - load + lt;
    end
end