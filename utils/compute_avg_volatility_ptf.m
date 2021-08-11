function [OptMVPtf, avg_vol] = compute_avg_volatility_ptf(sigma, mu)

    n = length(mu);
    N = 5000000;
    X = Sampling_simplex(n, N, 'RM');

    Vols = X' * sigma;
    Vols = sum(Vols' .* X);
    avg_vol = mean(Vols);
    
    A = [ eye(n) ; -eye(n)];
    b = [ ones(n,1) ; zeros(n,1)];
    Aeq = [ ones(1,n) ];
    beq = [ 1 ];
    
    EqPtf = (1 / n) * ones(n, 1);
    
    options = optimoptions('fmincon','MaxFunctionEvaluations', 20000, 'Display', 'off', 'ConstraintTolerance', 1e-10, 'MaxIterations',5000, 'OptimalityTolerance', 1e-10);
    
    %compute a feasible point || Ptf0 = (1 / NbAssets) * ones(NbAssets, 1);
    Ptf0 = fmincon(@(x) FindInitPtf(x, sigma, avg_vol), EqPtf, A, b, Aeq, beq, [], [], [], options);
    
    %compute the optimal mean-variance portfolio || OptMVPtf = GetOptMVPtf(mu, sigma, VolConst);
    OptMVPtf = fmincon(@(x) objOptMVPtf(x, mu), Ptf0, A, b, Aeq, beq, [], [], @(x) VarCons(x, sigma, avg_vol), options);



end