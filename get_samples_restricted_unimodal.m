function X = get_samples_restricted_unimodal(sigma, y0, a, N, R, r, psrf_target)
    
    if (nargin == 6)
        psrf_target = 1.02;
    end
    options = optimoptions('fmincon','MaxFunctionEvaluations', 20000, 'Display', 'off', 'ConstraintTolerance', 1e-10, 'MaxIterations',5000, 'OptimalityTolerance', 1e-10);

    n = length(y0);
    
    A = [ eye(n) ; -eye(n); R];
    b = [ ones(n,1) ; zeros(n,1); r];
    Aeq = [ ones(1,n) ];
    beq = [ 1 ];
    
    NN = null(Aeq);
    shift = (1/n)*ones(n,1);
    
    b = b - A * shift;
    A = A * NN;
    
    mu = y0 - shift;
    mu = (mu'*(sigma*NN))';
    q = 2;
    
    sigma = NN'*sigma*NN;
    
    sqrt_sum = sqrt(sum(A.^2,2));
    A = A ./ repmat(sqrt_sum, [1 n-1]);
    b = b ./ sqrt_sum;
    
    [xc,~] = get_cheb(A,b);
    %x0 = NN'*x0;
    
    x0 = fmincon(@(x) target_q_volatility_fun(x, sigma, mu, q), xc, A, b, [],[], [], [], [], options);
    %x0
    %NN*x0 + shift
    
    
    v_eq = xc - x0;
    v_eq = v_eq / norm(v_eq);
    
    l = (b-A*x0) ./ (A*v_eq);
    l = max(1./l);
    l = 1 / l;
    x0 = x0 + (0.1*l)*v_eq;
    
    %A*x0-b
    %sqrt_sum = sqrt(sum(A.^2,2))
    %A*x0-b
    [eps_step, x1, ~] = Initialize_hmc_leapfrog_Dual_Avg(A, b, x0, sigma, mu, a, q, 1000, 0.65);
    %A*x0-b
    %eps_step
    X_iter = hmc_leapfrog(A, b, x0, sigma, mu, a, q, 1500, eps_step/4);
    x0 = X_iter(:, 1500);
    n0 = 2000;
    N_total = 0;
    X = [];
    while(true)
        X_iter = hmc_leapfrog(A, b, x0, sigma, mu, a, q, n0, eps_step/4);
        N_total = N_total + n0;
        X = [X X_iter];
        psrf_iter = max(psrf(X'));
        if (psrf_iter <= psrf_target & N_total >= N)
            break;
        end
        x0 = X_iter(:, n0);
    end
    %X = hmc_leapfrog(A, b, x0, sigma, mu, a, q, N);
    
    X = NN * X + repmat(shift, [1 N_total]);

end


