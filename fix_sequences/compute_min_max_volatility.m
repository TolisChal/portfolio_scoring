function [min_vol, max_vol] = compute_min_max_volatility(sigma)

    n = size(sigma,1);
    EqPtf = (1/n)*ones(n,1);
    
    A = [ eye(n) ; -eye(n)];
    b = [ ones(n,1) ; zeros(n,1)];
    Aeq = [ ones(1,n) ];
    beq = [ 1 ];
    
    options = optimoptions('fmincon','MaxFunctionEvaluations', 20000, 'Display', 'off', 'ConstraintTolerance', 1e-10, 'MaxIterations',5000, 'OptimalityTolerance', 1e-10);
    
    Ptf = fmincon(@(x) volatility_fun(x, sigma), EqPtf, A, b, Aeq, beq, [], [], [], options);
    
    min_vol = Ptf'*sigma*Ptf;
    
    Ptf = fmincon(@(x) neg_volatility_fun(x, sigma), EqPtf, A, b, Aeq, beq, [], [], [], options);
   
    max_vol = Ptf'*sigma*Ptf;


end