function ptf = compute_erc_ptf(sigma)
    
    n = size(sigma,1);
    W=ones(n,1)/n;
    
    A = [ eye(n) ; -eye(n)];
    b = [ ones(n,1) ; zeros(n,1)];
    Aeq = [ ones(1,n) ];
    beq = [ 1 ];
    
    options = optimoptions('fmincon','MaxFunctionEvaluations', 20000, 'Display', 'off', 'ConstraintTolerance', 1e-10, 'MaxIterations',5000, 'OptimalityTolerance', 1e-10);
 
    f = @(W) var(W.*(sigma*W))*10^14; %Note: The 10^14 is there to increase accuracy
    ptf = fmincon(f,W,A, b, Aeq, beq, [], [], [], options);

end