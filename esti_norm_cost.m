function norm_const = esti_norm_cost(sigma, q, mu, N, alpha)

    n = length(mu);
    X = Sampling_simplex(n, N, 'RM');
    
    q_mu = q*(mu'*X);
    
    Sx = sigma*X;
    
    for i=1:N
   
        log_xSx(i) = X(:,i)' * Sx(:,i);
    
    end
    
    norm_const = mean(exp(- alpha * (log_xSx - q_mu)));

end