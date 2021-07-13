function ratios = compute_integral_ratios_fixed_q(sigma, mu, q, a_vals, as_dense, x0, N, R, r)
    
    n = length(mu);
    a_vals_sorted = sort(a_vals);
    as_dense_sorted = sort(as_dense);
    
    as_merged = sort([a_vals as_dense]);
    k = length(a_vals);
    ratios = zeros(1,k);
    
    a_vals_pos(1) = 1;
    for i = 1:k
        a_vals_pos(i+1) = find(as_merged == a_vals_sorted(i));
    end
    
    volume_ratio = Ali73(R, r);
    
    X1 = Sampling_simplex(n, N, 'RM');
    X2 = billiard_walk_low_dim(R, r, 2, N);
    
    a = as_dense_sorted(1);
    
    q_mu = q*(mu'*X1);
    Sx = sigma*X1;
    log_xSx = zeros(1,N);
    for j=1:N
        log_xSx(j) = X1(:,j)' * Sx(:,j);
    end
    const_int_num = mean(exp(-a*(log_xSx - q_mu)));
    
    q_mu = q*(mu'*X2);
    Sx = sigma*X2;
    log_xSx = zeros(1,N);
    for j=1:N
        log_xSx(j) = X2(:,j)' * Sx(:,j);
    end
    const_int_den = mean(exp(-a*(log_xSx - q_mu)));
    pos = a_vals_pos(1);
    for i = 1:k
        
        to_pos = a_vals_pos(i+1) - 1;
        a = a_vals_sorted(i);
        X = get_samples(sigma, mu, a, q, N, x0);
    
        sc = R * X;
        int_ratio = sum(sc < r) / size(X, 2)
        
        if (int_ratio >= 0.1)
            ratios(k+1-i) = int_ratio;
        else
            [int_ratio, const_int_num, const_int_den] = estimate_tele_integral_ratio(sigma, mu, q, as_merged, pos, to_pos, R, r, x0, const_int_num, const_int_den, volume_ratio, N);
            ratios(k+1-i) = int_ratio;
            pos = a_vals_pos(i+1);
        end
    end
end

