function values = compute_integrals_fixed_q(sigma, mu, q, a_vals, N, x0)
    
    n = length(mu);
    a_vals_sorted = sort(a_vals);
    
    X = Sampling_simplex(n, N, 'RM');
    a_prev = a_vals_sorted(1);
    k = length(a_vals_sorted);
    
    q_mu = q*(mu'*X);
    
    Sx = sigma*X;
    log_xSx = zeros(1,N);
    
    for i=1:N
        log_xSx(i) = X(:,i)' * Sx(:,i);
    end
    values = zeros(1,k);
    values(k) = mean(exp(-a_prev*(log_xSx - q_mu)));
    
    for i = 2:k
        a = a_vals_sorted(i);
        X = get_samples(sigma, mu, a_prev, q, N, x0);
        
        q_mu = q*(mu'*X);
    
        Sx = sigma*X;
        log_xSx = zeros(1,N);
    
        for j=1:N
            log_xSx(j) = X(:,j)' * Sx(:,j);
        end
    
        values(k+1-i) = mean(exp((-a + a_prev)*(log_xSx - q_mu)))* values(k+2-i);
        a_prev = a;
    end
        




end
