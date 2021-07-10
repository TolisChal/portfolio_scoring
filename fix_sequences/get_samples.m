function X = get_samples(sigma, mu, a, q, N, x0)
    
    n = length(mu);
    
    A = [ eye(n) ; -eye(n)];
    b = [ ones(n,1) ; zeros(n,1)];
    Aeq = [ ones(1,n) ];
    beq = [ 1 ];
    
    
    NN = null(Aeq);
    shift = (1/n)*ones(n,1);
    x0 = x0 - shift;
    
    b = b - A * shift;
    A = A * NN;
    
    sigma = NN'*sigma*NN;
    mu = (mu' * NN)';
    x0 = NN'*x0;
    
    [eps_step, x0, prob_step_in] = Initialize_hmc_leapfrog_Dual_Avg(A, b, X0, sigma, mu, a, q, N, 0.65);
    X = hmc_leapfrog(A, b, x0, sigma, mu, a, q, N, eps_step);
    %X = hmc_leapfrog(A, b, x0, sigma, mu, a, q, N);
    
    X = NN * X;

end


