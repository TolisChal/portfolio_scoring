function a_max = compute_max_a(sigma, mu, q, x0, N, a)
    
    n = length(mu);
    ratio = (1 + 1/(2*n));

    X = get_samples(sigma, mu, a, q, N, x0);
    
    dist_max_prev = max(sqrt(sum((X - repmat(x0, [1 size(X,2)])).^2,1)));
    iter = 1;
    while (iter <= 3)
        a = a*ratio;
        X = get_samples(sigma, mu, a, q, N, x0);
        
        dist_max = max(sqrt(sum((X - repmat(x0, [1 size(X,2)])).^2,1)));
        %dist_max_prev
        %abs(dist_max - dist_max_prev)/dist_max_prev
        if (abs(dist_max - dist_max_prev)/dist_max_prev < 0.05)
            a_max = a;
            return
        end
        dist_max_prev = dist_max;
        
        iter = iter + 1;
    end
    a_max = a;

end