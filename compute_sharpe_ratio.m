function sh_ratio = compute_sharpe_ratio(R, x0, sigma)
    
    n = length(x0);
    EqPtf = ones(n,1)/n;
    sh_ratio = (R*x0 - R*EqPtf) / sqrt(x0'*sigma*x0);

end