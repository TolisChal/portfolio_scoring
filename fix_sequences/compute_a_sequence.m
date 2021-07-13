function [as, as_dense, samples] = compute_a_sequence(sigma, mu, q, max_vol, M, x0, N)

    n = length(mu);
    
    a_min = (1 / (max_vol + q*max(mu)));    
    a_max = 7*n*n;
    
    a_max = compute_max_a(sigma, mu, q, x0, 50000, a_max);
    
    [as_dense, samples] = compute_a_dense_sequence(sigma, mu, q, a_max, a_min, x0, N);
    as_dense
    %X = get_samples(sigma, mu, a_min, q, N, x0);
    %get_samples(sigma, mu, a, q, 5000, x0);
    
    dist_total = 0;
    k = length(as_dense);
    dists_dense(1) = 0;
    for i = 1 : 1:(k-1)
        dists_dense(i+1) = estimate_L2_norm(sigma, mu, q, as_dense(i), as_dense(i+1), samples{i})
        dist_total = dist_total + dists_dense(i);
    end

    %cum_dists_dense = cumsum(dists_dense);
    step = dist_total / (M+1);
    %dist_total;
    
    dists = 0:step:dist_total;
    dists = dists(2:end-1)
    
    a = as_dense(1);
    samples = {};
    
    %X = samples{1};
    X = get_samples(sigma, mu, a, q, N, x0);
    counter = 1;
    as_dense
    dist_prev = 0;
    for i = 1:length(dists)
        %indx = find(cum_dists_dense > i);
        %indx = indx(1);
        d = dists(i) - dist_prev;
        a_next = compute_next_a(sigma, mu, a, as_dense(end), q, X, d);
        dist_prev = dists(i);
        
        X = get_samples(sigma, mu, a_next, q, N, x0);
        samples{i} = X;
        
        as(counter) = a_next;
        counter = counter + 1;
        
        a = a_next;
    end
end