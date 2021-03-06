function a_min = compute_min_a(sigma, mu, q, N, a0_min, x0)

    n = length(mu);
    X = Sampling_simplex(n, N, 'RM');
    
    q_mu = q*(mu'*X);
    
    Sx = sigma*X;
    log_xSx = zeros(1,N);
    
    for i=1:N
        log_xSx(i) = X(:,i)' * Sx(:,i);
    end
    
    ratio = 1 - (1/n);
    
    a_min = a0_min;
    
    while(true)
        
        a_min = a_min*ratio;
        
        max_val = exp(-a_min*((x0'*sigma*x0)-q*(mu'*x0)));
        
        vals = exp(-a_min*(log_xSx - q_mu));
        
        min_val = min(vals);
        %abs(max_val - min_val) / abs(max_val)
        if (abs(max_val / min_val) < 100) 
            return
        end
        
    end



end