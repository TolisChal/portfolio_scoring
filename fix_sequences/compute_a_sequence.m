function as = compute_a_sequence(sigma, mu, q, M, max_vol)

    n = length(mu);

    a_min = 1 / (max_vol + q*max(mu));    
    a_max = 4*n;
    
    X = get_samples(sigma, mu, a_min, q, N);
    
    dist_max = compute_L2_norm(sigma, mu, a_min. a_max, q, X);
    dist_min = 0;

    step = (dist_min - dist_max) / (M+1);
    
    dists = dist_min:step:dist_max;
    dists = dists(2:end-1);
    
    a = a_min;
    N = 10000;
    for i = dists
    
        X = get_samples(sigma, mu, a, q, N);
        
        a_next = compute_next_a(sigma, mu, a, q, X, i);
    
    end
end