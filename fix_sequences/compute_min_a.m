function a_min = compute_min_a(sigma, mu, q, N, a0_min)

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
        
        vals = exp(-a_min*(log_xSx - q_mu));
        
        abs(max(vals) - min(vals)) / abs(max(vals))
        if (abs(max(vals) - min(vals)) / abs(max(vals)) < 0.5) 
            return
        end
        
    end



end