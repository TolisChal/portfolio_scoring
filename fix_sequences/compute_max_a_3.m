function a_max = compute_max_a_3(sigma, mu, q, N, a0_max, x0)

    n = length(mu);
    ratio = 1 + (1/n);
    
    a0_max = a0_max*(1 - (1/n));
    a_max = a0_max;
    
    curr_fn = 2;
    
    while(curr_fn < 3.35*(10^8))
        
        a_prev = a_max;
        a_max = a_max*ratio;
        
        X = get_samples(sigma, mu, a_max, q, N, x0);
    
        q_mu = q*(mu'*X);
    
        Sx = sigma*X;
        log_xSx = zeros(1,size(X,2));
    
        for i=1:size(X,2)
            log_xSx(i) = X(:,i)' * Sx(:,i);
        end
        
        curr_fn = mean(exp((- a_max + a_prev) * (log_xSx - q_mu)));
        
    end
    
end