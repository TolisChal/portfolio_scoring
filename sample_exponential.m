function X = sample_exponential(c, T, N, A_add, b_add)

    n = length(c);
    
    psrf_target = 1.02;
    
    A = [ eye(n) ; -eye(n)];
    b = [ ones(n,1) ; zeros(n,1)];
    Aeq = [ ones(1,n) ];
    beq = [ 1 ];

    if (nargin > 3)
        A = [A; A_add];
        b = [b; b_add];
    end
    
    NN = null(Aeq);
    shift = (1/n)*ones(n,1);
    %x0 = x0 - shift;
    
    b = b - A * shift;
    A = A * NN;
    
    c = (c' * NN)';
    
    [x0, ~]=get_cheb(A,b);
    
    sqrt_sum = sqrt(sum(A.^2,2));
    A = A ./ repmat(sqrt_sum, [1 n-1]);
    b = b ./ sqrt_sum;
    %sqrt_sum = sqrt(sum(A.^2,2))
    %A*x0-b
    [eps_step, x0, ~] = Initialize_hmc_exp_leapfrog_Dual_Avg(A, b, x0, c, T, 500, 0.65);
    %eps_step
    %r
    n0 = 2000;
    N_total = 0;
    X = [];
    while(true)
        X_iter = hmc_exp_leapfrog(A, b, x0, c, T, n0, eps_step);
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