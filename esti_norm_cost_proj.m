function norm_const = esti_norm_cost_proj(sigma, q, mu, N, alpha)

    n = length(mu);
    
    Aeq = ones(1,n);
    
    NN = null(Aeq);
    
    prod(svd(NN))
    
    A = [-eye(n); ones(1,n)];
    b = [zeros(n,1); 1];
    
    shift = (1/n)*ones(n,1);
    
    b = b - A * shift;
    A = A * NN;
    
    x = zeros(n-1,1);
    X = billiard_walk(A, b, x, 3, N, 2);
    X = NN*X + repmat(shift, [1 N]);
    sum(X)
    size(X)
    
    %sigma = NN' * sigma * NN;
    %mu = (mu'*NN)';
    size(sigma)
    size(mu)
    
    %X = Sampling_simplex(n, N, 'RM');
    
    q_mu = q*(mu'*X);
    
    Sx = sigma*X;
    
    for i=1:N
   
        log_xSx(i) = X(:,i)' * Sx(:,i);
    
    end
    
    norm_const = mean(exp(- alpha * (log_xSx - q_mu)));

end