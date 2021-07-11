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
    
    size(A)
    sqrt_sum = sqrt(sum(A.^2,2));
    A = A ./ repmat(sqrt_sum, [1 n-1]);
    b = b ./ sqrt_sum;
    %sqrt_sum = sqrt(sum(A.^2,2))
    
    [eps_step, x0, ~] = Initialize_hmc_leapfrog_Dual_Avg(A, b, x0, sigma, mu, a, q, 1000, 0.65);
    eps_step
    X = hmc_leapfrog(A, b, x0, sigma, mu, a, q, N, eps_step/4);
    %X = hmc_leapfrog(A, b, x0, sigma, mu, a, q, N);
    
    X = NN * X + repmat(shift, [1 N]);

end


