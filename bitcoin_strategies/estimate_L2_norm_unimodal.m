function dist = estimate_L2_norm_unimodal(sigma, y0, a1, a2, X)

    n = length(y0);
    N = size(X, 2);
    
    X = X - repmat(y0, [1, N]);
    Sx = sigma*X;
    log_xSx = zeros(1,N);
    
    for i=1:N
   
        log_xSx(i) = X(:,i)' * Sx(:,i);
    
    end

    dist = mean(exp((- a1 + a2) * (log_xSx))) * mean(exp((a1 - a2) * (log_xSx)));
    
end

