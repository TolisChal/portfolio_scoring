function ratios = compute_integral_ratios_fixed_q_all(sigma, mu, q, a_vals, x0, N, R, r)
    
    n = length(mu);
    a_vals_sorted = sort(a_vals);

    k = length(a_vals);
    ratios = zeros(1,k);
    
    volume_ratio = Ali73(R, r);
    
    X1 = Sampling_simplex(n, N, 'RM');
    X2 = billiard_walk_low_dim(R, r, 2, N);
    
    a = a_vals_sorted(1);
    
    q_mu = q*(mu'*X1);
    Sx = sigma*X1;
    log_xSx = zeros(1,N);
    for j=1:N
        log_xSx(j) = X1(:,j)' * Sx(:,j);
    end
    const_int_den = mean(exp(-a*(log_xSx - q_mu)));
    
    q_mu = q*(mu'*X2);
    Sx = sigma*X2;
    log_xSx = zeros(1,N);
    for j=1:N
        log_xSx(j) = X2(:,j)' * Sx(:,j);
    end
    const_int_num = mean(exp(-a*(log_xSx - q_mu)));
    pos = 1;
    
    ratios(k) = volume_ratio * (const_int_num / const_int_den);
    
    for i = 2:(k-1)
        
        to_pos = i-1;
        a = a_vals_sorted(i);
        X = get_samples(sigma, mu, a, q, N, x0);
    
        sc = R * X;
        int_ratio = sum(sc < r) / size(X, 2)
        
        if (int_ratio >= 0.2)
            ratios(k+1-i) = int_ratio;
        else
            [int_ratio, const_int_num, const_int_den] = estimate_tele_integral_ratio(sigma, mu, q, a_vals_sorted, pos, to_pos, R, r, x0, const_int_num, const_int_den, volume_ratio, N);
            int_ratio
            ratios(k+1-i) = int_ratio;
            pos = i;
        end
    end
end

