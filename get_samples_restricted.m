function X = get_samples_restricted(sigma, mu, a, q, N, R, r)
    
    n = length(mu);
    
    A = [ eye(n) ; -eye(n); R];
    b = [ ones(n,1) ; zeros(n,1); r];
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
    
    sigma = NN'*sigma*NN;
    mu = (mu' * NN)';
    x0 = NN'*x0;
    
    
    sqrt_sum = sqrt(sum(A.^2,2));
    A = A ./ repmat(sqrt_sum, [1 n-1]);
    b = b ./ sqrt_sum;
    %sqrt_sum = sqrt(sum(A.^2,2))
    
    [eps_step, x0, ~] = Initialize_hmc_leapfrog_Dual_Avg(A, b, x0, sigma, mu, a, q, 1000, 0.65);
    %eps_step
    X = hmc_leapfrog(A, b, x0, sigma, mu, a, q, N, eps_step/4);
    %X = hmc_leapfrog(A, b, x0, sigma, mu, a, q, N);
    
    X = NN * X + repmat(shift, [1 N]);

end


