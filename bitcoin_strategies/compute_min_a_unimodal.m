function a_min = compute_min_a_unimodal(sigma, y0, N, a0_min)

    n = length(y0);
    X = Sampling_simplex(n, N, 'RM');
    
    X = X - repmat(y0, [1, N]);
    Sx = sigma*X;
    log_xSx = zeros(1,N);
    
    for i=1:N
        log_xSx(i) = X(:,i)' * Sx(:,i);
    end
    
    ratio = 1 - (1/n);
    
    a_min = a0_min;
    
    while(true)
        
        a_min = a_min*ratio;
        
        max_val = 1;
        
        vals = exp(-a_min*(log_xSx));
        
        min_val = min(vals);
        %abs(max_val - min_val) / abs(max_val)
        if (abs(max_val / min_val) < 100) 
            return
        end
        
    end



end