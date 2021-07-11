function a_next = compute_next_a(sigma, mu, a_min, a_max, q, X, dist)

    a1_iter = a_min;
    a2_iter = a_max;
    
    while(true)
        
        a_med = (a1_iter + a2_iter);
        
        d = estimate_L2_norm(sigma, mu, q, a1_iter, a_med, X);
        
        if (abs(d-dist)/dist < 0.01)
            a_next = a_med;
            return
        end
        
    end

end

