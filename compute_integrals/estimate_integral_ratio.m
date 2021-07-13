function int_ratio = estimate_integral_ratio(sigma, mu, a, q, N, x0, R, r)

    X = get_samples(sigma, mu, a, q, N, x0);
    
    sc = R * X;
    int_ratio = sum(sc < r) / N;
    
    if (int_ratio > 0.05)
        return
    end

end