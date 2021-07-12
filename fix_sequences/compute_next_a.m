function a_next = compute_next_a(sigma, mu, a, a_min, q, X, dist)

    a1_iter = a_min;
    a2_iter = a;
    
    iter = 1;
    while(iter < 100)
        
        a_med = (a1_iter + a2_iter)/2;
        
        d = estimate_L2_norm(sigma, mu, q, a, a_med, X);
        %dist
        
        if (abs(d-dist)/dist < 0.02)
            a_next = a_med;
            return
        elseif (d < dist)
            a2_iter = a_med;
        else
            a1_iter = a_med;
        end
        
        iter = iter + 1;
    end
    a_next = a_med;
    d
    dist
end

