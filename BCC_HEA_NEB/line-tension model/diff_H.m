function f=diff_H(para, tau_xy, pb,z0,position)
    % para(1) -> alpha
    % para(2) -> beta
    % f(1) -> \partial H / \partial alpha
    % f(2) -> \partial H / \partial beta
    constants;
    tau_xy = tau_xy / converter; % convert MPa to ev/A^3
    alpha = para(1);
    beta = para(2);
    f(1) = alpha;
    f(2) = alpha;

    X = (-N):1:N; % position along dislocation line
    %z0 = 0;%1/2/pi*(asin(tau_xy*b^2*Lp/DeltaU/pi));
    Z = zeros(1,2*N+1);
    for n = 1:(2*N+1)
        Z(n) = (0.5 * tanh(alpha * (X(n)-position + beta)) - 0.5 * tanh(alpha * (X(n)-position - beta)))*(1-z0)+z0;
    end

    for n = 1:(2*N+1)
        if n ~= 2*N+1
            subs_p = Z(n) - Z(n+1);
        else
            subs_p = Z(n) - Z(1);
        end
        if n ~= 1
            subs_m = Z(n) - Z(n-1);
        else
            subs_m = Z(n) - Z(2*N+1);
        end
        hpartialz = pb(n)*pi*sin(2*pi*Z(n)/period)/period - tau_xy*b^2*Lp + Gamma*Lp^2/b * subs_p +  Gamma* Lp^2/b * subs_m;
        zpartiala = 0.5 * ( sech(alpha * (X(n)-position + beta)) * sech(alpha * (X(n)-position + beta)) * (X(n)-position + beta) - sech(alpha * (X(n)-position - beta)) * sech(alpha * (X(n)-position - beta)) * (X(n)-position - beta) );
        zpartialb = 0.5 * ( sech(alpha * (X(n)-position + beta)) * sech(alpha * (X(n)-position + beta)) * alpha + sech(alpha * (X(n)-position - beta)) * sech(alpha * (X(n)-position - beta)) * alpha );
        f(1) = f(1) + hpartialz * zpartiala * (1-z0);
        f(2) = f(2) + hpartialz * zpartialb * (1-z0);
    end
    f(1)=f(1)-alpha;
    f(2)=f(2)-alpha;
    
end