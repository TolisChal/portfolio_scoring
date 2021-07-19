function a_max = compute_max_a_2_unimodal(sigma, y0, N, a0_max)

    n = length(y0);
    X = Sampling_simplex(n, N, 'RM');
    
    X = X - repmat(y0, [1, N]);
    Sx = sigma*X;
    log_xSx = zeros(1,N);
    
    for i=1:N
        log_xSx(i) = X(:,i)' * Sx(:,i);
    end
    
    ratio = 1 + (1/n);
    
    a_max = a0_max;
    
    while(true)
        
        a_max = a_max*ratio;
        
        max_val = 1;
        
        vals = exp(-a_max*(log_xSx));
        
        min_val = min(vals);
        abs(max_val / min_val)
        if (abs(max_val / min_val) > 10^100) 
            return
        end
        
    end



end