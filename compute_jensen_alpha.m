function Ja = compute_jensen_alpha(returns, R, x0, sigma, mu)

    n = length(x0);
    EqPtf = ones(n,1)/n;

    Ret = R*x0;
    
    Ret_eq = (returns*EqPtf);
    ret_x0 = (returns*x0);
    beta = cov(ret_x0,Ret_eq);
    beta = beta(1,2) / var(Ret_eq);
    
    market_return = R*EqPtf;
    
    A = [ eye(n) ; -eye(n)];
    b = [ ones(n,1) ; zeros(n,1)];
    Aeq = [ ones(1,n) ];
    beq = [ 1 ];
    
    options = optimoptions('fmincon','MaxFunctionEvaluations', 20000, 'Display', 'off', 'ConstraintTolerance', 1e-10, 'MaxIterations',5000, 'OptimalityTolerance', 1e-10);
    
    Ptf = fmincon(@(x) volatility_fun(x, sigma), EqPtf, A, b, Aeq, beq, [], [], [], options);
    risk_free_rate = mu'*Ptf;
    
    Ja = Ret - (risk_free_rate + beta *(market_return - risk_free_rate));
    
end