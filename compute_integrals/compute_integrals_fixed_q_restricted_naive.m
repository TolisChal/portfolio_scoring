function values = compute_integrals_fixed_q_restricted_naive(sigma, mu, q, a_vals, N, R, r)
    
    n = length(mu);
    a_vals_sorted = sort(a_vals);
    
    X = billiard_walk_low_dim(R, r, 2, N);
    
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