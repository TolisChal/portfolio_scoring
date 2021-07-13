function X = get_samples_restricted(sigma, mu, a, q, N, R, r, psrf_target)
    
    if (nargin == 7)
        psrf_target = 1.02;
    end

    n = length(mu);
    
    A = [ eye(n) ; -eye(n); R];
    b = [ ones(n,1) ; zeros(n,1); r];
    Aeq = [ ones(1,n) ];
    beq = [ 1 ];

    NN = null(Aeq);
    shift = (1/n)*ones(n,1);

    b = b - A * shift;
    A = A * NN;
    
    mu = (mu - (2/q)*(shift'*sigma)');
    
    sigma = NN'*sigma*NN;
    mu = (mu' * NN)';

    [x0,~] = get_cheb(A,b);
    
    %A*x0-b
    
    sqrt_sum = sqrt(sum(A.^2,2));
    A = A ./ repmat(sqrt_sum, [1 n-1]);
    b = b ./ sqrt_sum;
    %sqrt_sum = sqrt(sum(A.^2,2))
    
    [eps_step, x0, ~] = Initialize_hmc_leapfrog_Dual_Avg(A, b, x0, sigma, mu, a, q, 1000, 0.65);
    %sum((A*x0-b)>0)
    %eps_step
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
    %X = hmc_leapfrog(A, b, x0, sigma, mu, a, q, N, eps_step/4);
    %X = hmc_leapfrog(A, b, x0, sigma, mu, a, q, N);
    
    %n_out = 0;
    %for i=1:N
    %    if (sum((A*X(:,i)-b)>0) > 0)
    %        n_out = n_out + 1;
    %    end
    %end
    %n_out
    
    X = NN * X + repmat(shift, [1 N_total]);

end


