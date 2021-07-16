function T_max = compute_max_T(c, T_min)

    n = length(c);
    X = Sampling_simplex(n, N, 'RM');
    
    log_vals = -c'*X;
    
    ratio = 1 + (1/n);
    
    T_max = T_min;
    
    while(true)
        
        T_max = T_max*ratio;
        
        vals = exp(log_vals / T_max);
        
        abs(max(vals) - min(vals)) / abs(max(vals))
        if (abs(max(vals) - min(vals)) / abs(max(vals)) < 0.1) 
            return
        end
        
    end


end