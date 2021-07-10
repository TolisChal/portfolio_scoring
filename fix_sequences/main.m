clear all; close all; clc;

%Read the data
T = readtable('asset_returns_10.csv');
AssRets = table2array(T) / 100;

IR = readtable('IR.csv');
IR = table2array(IR) / 100;

% number of assets and number asset returns
[N, NbAssets] = size(AssRets);

EqPtf = (1 / NbAssets) * ones(NbAssets, 1); %Equally weighted portfolio
RndPtfs = Sampling(NbAssets, 10^6, 'RK'); %uniformly distributed portfolios
W = 120; %sliding window

%shift the assets' returns by the return of the equally weighted portfolio
%EqPtf_rets = AssRets * EqPtf;
%AssRets = AssRets - repmat(EqPtf_rets, [1 NbAssets]);
AssRets = AssRets - repmat(IR(:,4), [1 NbAssets]);

%The linear constriants of the canonical simplex (the set of all
%portfolios)
A = [ eye(NbAssets) ; -eye(NbAssets)];
b = [ ones(NbAssets,1) ; zeros(NbAssets,1)];
Aeq = [ ones(1,NbAssets) ];
beq = [ 1 ];
c
options2 = optimoptions('fmincon','MaxFunctionEvaluations', 20000, 'Display', 'off', 'ConstraintTolerance', 1e-10, 'MaxIterations',5000, 'OptimalityTolerance', 1e-10);

cell_of_EF_figures = {};

%for i=1:1
for i=1:(N-W)
    
    disp(' ')
    disp(strcat('iteration = ', num2str(i)))
    
    tic
    
    AssRets_W = AssRets(i:(i+W-1), :); %asset returns of this sliding window
    Log_AssRets_W = log(AssRets_W + 1);
    
    %estimate the covariance and the mean
    [sigma, shr] = covCor(Log_AssRets_W);
    mu = mean(Log_AssRets_W)';
    
    %generate normally distributed asset returns
    %RndAssRets = mvnrnd(mu, sigma, 10^5);
    
    %compute the average volatility of the set of random portfolios
    Vols = RndPtfs' * sigma;
    Vols = sum(Vols' .* RndPtfs);
    VolConst = mean(Vols);
    
    %compute a feasible point || Ptf0 = (1 / NbAssets) * ones(NbAssets, 1);
    Ptf0 = fmincon(@(x) FindInitPtf(x, sigma, VolConst), EqPtf, A, b, Aeq, beq, [], [], [], options);
    
    %compute the optimal mean-variance portfolio || OptMVPtf = GetOptMVPtf(mu, sigma, VolConst);
    OptMVPtf = fmincon(@(x) objOptMVPtf(x, mu), Ptf0, A, b, Aeq, beq, [], [], @(x) VarCons(x, sigma, VolConst), options);
    disp('OptMVPtf computed')
    
    %compute the optimal Perf_A portfolio || OptScPtf_A = GetOptScPtf_A(sigma, VolConst, RndAssRets, NbAssets);
    OptScPtf_A = fmincon(@(x) objOptScPtf_A(x, AssRets_W),Ptf0, A, b, Aeq, beq, [], [], @(x) VarCons(x, sigma, VolConst), options2);
    disp('OptScPtf_A computed')
    
    %compute the optimal Perf_B portfolio || OptScPtf_B = GetOptScPtf_B(sigma, VolConst, RndAssRets, NbAssets);
    OptScPtf_B = fmincon(@(x) objOptScPtf_B(x, AssRets_W),Ptf0, A, b, Aeq, beq, [], [], @(x) VarCons(x, sigma, VolConst), options2);
    disp('OptScPtf_B computed')
    
    %compute the optimal Perf_C portfolio || OptScPtf_C = GetOptScPtf_C(sigma, VolConst, RndAssRets, NbAssets);
    OptScPtf_C = fmincon(@(x) objOptScPtf_C(x, AssRets_W),Ptf0, A, b, Aeq, beq, [], [], @(x) VarCons(x, sigma, VolConst), options2);
    disp('OptScPtf_C computed')
    
    %compute the optimal Perf_D portfolio || OptScPtf_D = GetOptScPtf_D(sigma, VolConst, RndAssRets, NbAssets);
    OptScPtf_D = fmincon(@(x) objOptScPtf_D(x, AssRets_W),Ptf0, A, b, Aeq, beq, [], [], @(x) VarCons(x, sigma, VolConst), options2);
    disp('OptScPtf_D computed')
    
    AssRets_next = AssRets(i + W, :);
    
    OptMVPtf_scores(i) = Ali73_vecRet(AssRets_next, OptMVPtf);
    OptScPtf_A_scores(i) = Ali73_vecRet(AssRets_next, OptScPtf_A);
    OptScPtf_B_scores(i) = Ali73_vecRet(AssRets_next, OptScPtf_B);
    OptScPtf_C_scores(i) = Ali73_vecRet(AssRets_next, OptScPtf_C);
    OptScPtf_D_scores(i) = Ali73_vecRet(AssRets_next, OptScPtf_D);
    
    OptMVPtf_returns(i) = AssRets_next * OptMVPtf;
    OptScPtf_A_returns(i) = AssRets_next * OptScPtf_A;
    OptScPtf_B_returns(i) = AssRets_next * OptScPtf_B;
    OptScPtf_C_returns(i) = AssRets_next * OptScPtf_C;
    OptScPtf_D_returns(i) = AssRets_next * OptScPtf_D;
    
    Ptfs_Opt{i} = OptMVPtf;
    Ptfs_A{i} = OptScPtf_A;
    Ptfs_B{i} = OptScPtf_B;
    Ptfs_C{i} = OptScPtf_C;
    Ptfs_D{i} = OptScPtf_D;
    
    means{i} = mu;
    cov_mats{i} = sigma;
    
    %Plot_Efficient_Frontier(mu, sigma, OptMVPtf, OptScPtf_A, OptScPtf_B, OptScPtf_C, OptScPtf_D, i);
   
    iter_runtime = toc;
    runtimes(i) = iter_runtime;
    
    if (mod(i,50)==0 || i==(N-W))
        save('OptMVPtf_scores.mat', 'OptMVPtf_scores') 
        save('OptScPtf_A_scores.mat', 'OptScPtf_A_scores')
        save('OptScPtf_B_scores.mat', 'OptScPtf_B_scores')
        save('OptScPtf_C_scores.mat', 'OptScPtf_C_scores')
        save('OptScPtf_D_scores.mat', 'OptScPtf_D_scores')
    
        save('OptMVPtf_returns.mat', 'OptMVPtf_returns')
        save('OptScPtf_A_returns.mat', 'OptScPtf_A_returns')
        save('OptScPtf_B_returns.mat', 'OptScPtf_B_returns')
        save('OptScPtf_C_returns.mat', 'OptScPtf_C_returns')
        save('OptScPtf_D_returns.mat', 'OptScPtf_D_returns')
    
        save('Ptfs_Opt.mat', 'Ptfs_Opt')
        save('Ptfs_A.mat', 'Ptfs_A')
        save('Ptfs_B.mat', 'Ptfs_B')
        save('Ptfs_C.mat', 'Ptfs_C')
        save('Ptfs_D.mat', 'Ptfs_D')
    
        save('means.mat', 'means')
        save('cov_mats.mat', 'cov_mats')
    
        save('runtimes.mat', 'runtimes')
    end
    
    disp(strcat('iteration completed, total runtime = ', num2str(iter_runtime), ' seconds'))
    disp(' ')
    
end