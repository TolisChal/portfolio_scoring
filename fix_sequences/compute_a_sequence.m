function as = compute_a_sequence(sigma, mu, q, M, max_vol, x0, N)

    n = length(mu);

    a_min = (1 / (max_vol + q*max(mu)))*2;    
    a_max = 40*n;
    
    X = get_samples(sigma, mu, a_min, q, N, x0);
    %get_samples(sigma, mu, a, q, 5000, x0);
    
    dist_max = estimate_L2_norm(sigma, mu, q, a_min, a_max, X);
    
    dist_min = 0;

    step = (dist_max - dist_min) / (M+1);
    
    dists = dist_min:step:dist_max;
    dists = dists(2:end-1);
    
    a = a_min;
    N = 10000;
    counter = 1;
    for i = dists

        a_next = compute_next_a(sigma, mu, a, a_max, q, X, i);
        
        X = get_samples(sigma, mu, a_next, q, N, x0);
        
        as(counter) = a_next;
        counter = counter + 1;
        
        a = a_next;
    end
end