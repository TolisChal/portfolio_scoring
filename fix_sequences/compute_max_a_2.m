function a_max = compute_max_a_2(sigma, mu, q, N, a0_max, x0)

    n = length(mu);
    X = Sampling_simplex(n, N, 'RM');
    
    q_mu = q*(mu'*X);
    
    Sx = sigma*X;
    log_xSx = zeros(1,N);
    
    for i=1:N
        log_xSx(i) = X(:,i)' * Sx(:,i);
    end
    
    ratio = 1 + (1/n);
    
    a_max = a0_max;
    
    while(true)
        
        a_max = a_max*ratio;
        
        max_val = exp(-a_max*((x0'*sigma*x0)-q*(mu'*x0)));
        
        vals = exp(-a_max*(log_xSx - q_mu));
        
        min_val = min(vals);
        %abs(max_val / min_val)
        if (abs(max_val / min_val) > 10^100) 
            return
        end
        
    end



end