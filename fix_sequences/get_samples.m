function X = get_samples(sigma, mu, a, q, N, x0, psrf_target)
    
    if (nargin == 6)
        psrf_target = 1.02;
    end

    n = length(mu);
    
    A = [ eye(n) ; -eye(n)];
    b = [ ones(n,1) ; zeros(n,1)];
    Aeq = [ ones(1,n) ];
    beq = [ 1 ];
    
    x_eq = (1/n)*ones(n,1);
    v_eq = x_eq - x0;
    v_eq = v_eq / norm(v_eq);
    
    l = (b-A*x0) ./ (A*v_eq);
    l = max(1./l);
    l = 1 / l;
    x0 = x0 + (0.1*l)*v_eq;
    %sum(x0)
    
    NN = null(Aeq);
    shift = (1/n)*ones(n,1);
    x0 = x0 - shift;
    
    b = b - A * shift;
    A = A * NN;
    
    mu = (mu - (2/q)*(shift'*sigma)');
    
    sigma = NN'*sigma*NN;
    mu = (mu' * NN)';
    x0 = NN'*x0;
    
    
    sqrt_sum = sqrt(sum(A.^2,2));
    A = A ./ repmat(sqrt_sum, [1 n-1]);
    b = b ./ sqrt_sum;
    %sqrt_sum = sqrt(sum(A.^2,2))
    %A*x0-b
    [eps_step, x1, ~] = Initialize_hmc_leapfrog_Dual_Avg(A, b, x0, sigma, mu, a, q, 1000, 0.65);
    %A*x0-b
    %eps_step
    X_iter = hmc_leapfrog(A, b, x0, sigma, mu, a, q, 1000, eps_step/4);
    x0 = X_iter(:, 1000);
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
        if (N_total >= 100000)
            break;
        end
        x0 = X_iter(:, n0);
    end
    %X = hmc_leapfrog(A, b, x0, sigma, mu, a, q, N);
    
    X = NN * X + repmat(shift, [1 N_total]);

end


