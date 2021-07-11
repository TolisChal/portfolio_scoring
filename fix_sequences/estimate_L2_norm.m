function dist = estimate_L2_norm(sigma, mu, q, a1, a2, X)

    n = length(mu);
    N = size(X, 2);
    
    q_mu = q*(mu'*X);
    
    Sx = sigma*X;
    log_xSx = zeros(1,N);
    
    for i=1:N
   
        log_xSx(i) = X(:,i)' * Sx(:,i);
    
    end

    dist = mean(exp((- a1 + a2) * (log_xSx - q_mu))) * mean(exp((a1 - a2) * (log_xSx - q_mu)));
    
end

