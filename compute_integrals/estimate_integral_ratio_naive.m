function int_ratio = estimate_integral_ratio_naive(sigma, mu, a, q, N, R, r)
    
    n = length(mu);
    
    X = billiard_walk_low_dim(R, r, 2, N);
 
    q_mu = q*(mu'*X);
    
    Sx = sigma*X;
    log_xSx = zeros(1,N);
    
    for i=1:N
        log_xSx(i) = X(:,i)' * Sx(:,i);
    end
    
    v1 = mean(exp(-a*(log_xSx - q_mu)));
    
    X = Sampling_simplex(n, N, 'RM');
 
    q_mu = q*(mu'*X);
    
    Sx = sigma*X;
    log_xSx = zeros(1,N);
    
    for i=1:N
        log_xSx(i) = X(:,i)' * Sx(:,i);
    end
    
    v2 = mean(exp(-a*(log_xSx - q_mu)));

    int_ratio = (v1 / v2)*Ali73(R,r);
    
end


