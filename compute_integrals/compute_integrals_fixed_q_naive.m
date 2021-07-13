function values = compute_integrals_fixed_q_naive(sigma, mu, q, a_vals, N)
    
    n = length(mu);
    a_vals_sorted = sort(a_vals);
    
    X = Sampling_simplex(n, 100000, 'RM');
    
    k = length(a_vals_sorted);
    values = zeros(1,k);
    
    for j = 1:k
        a = a_vals_sorted(j);
        
        q_mu = q*(mu'*X);
    
        Sx = sigma*X;
        log_xSx = zeros(1,100000);
    
        for i=1:100000
            log_xSx(i) = X(:,i)' * Sx(:,i);
        end
    
        values(k+1-j) = mean(exp(-a*(log_xSx - q_mu)));
    end

end
