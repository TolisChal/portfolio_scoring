function T_max = compute_a_max_exponential(w, N, T0_max)
    
    n = length(w);
    
    T_max = T0_max;
    Eqptf = ones(n,1) / n;
    
    ratio = 1 + (1/n);
    xc = ones(n,1)/n;
    
    while(true)
        
        T_max = T_max * ratio;
        X = sample_exponential(-w, T_max, N, xc); %, A_add, b_add)
        xc = mean(X,2);
        
        sum(abs(xc - Eqptf))
        if (sum(abs(xc - Eqptf)) < 0.03)
            return
        end
        
    end

end