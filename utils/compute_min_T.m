function T_min = compute_min_T(c, T_max)

    n = length(c);
    X = Sampling_simplex(n, N, 'RM');
    
    log_vals = -c'*X;
    
    ratio = 1 - (1/n);
    
    T_min = T_max;
    
    while(true)
        
        T_min = T_min*ratio;
        
        vals = exp(log_vals / T_min);
        
        abs(max(vals) - min(vals)) / abs(max(vals))
        if (abs(max(vals) - min(vals)) / abs(max(vals)) > 100) 
            return
        end
        
    end


end