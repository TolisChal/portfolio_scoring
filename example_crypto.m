addpath('./fix_sequences')
addpath('./mcmc_diagnostics')
addpath('./compute_integrals')

close all

load('crypto_12_historical.mat')

days = 400:460;
Rets = Returns_12_coins(days, :);

d = 12;
M = 3;

N = 10000;

mu = (prod(1 + Rets) - 1)';

[sigma, ~] = covCor(Rets);

[min_vol, max_vol] = compute_min_max_volatility(sigma);
[qs, ptfs, q_min, q_max, vols] = compute_q_sequence(sigma, mu, min_vol, max_vol, 3);

%q=qs(2);
%a=2;

%X = Sampling_simplex(d, 10000, 'RM');

%vals = Evaluate_density(X, sigma, mu, a, q);

%scatter3(X(1,:),X(2,:),X(3,:),25,vals,'.')


%Ptf = Sampling_simplex(d,1, 'RM');

%R = randn(1, d);
%r = R * Ptf;

%[OptMVPtf, avg_vol] = compute_avg_volatility_ptf(sigma, mu);
[OptMVPtf, avg_vol] = compute_2avg_volatility_ptf(sigma, mu);

R = Returns_12_coins(days(end)+1, :);
r = R * OptMVPtf;

integral_ratios = {};
iter = 1;
%as_all = zeros(5, 5);
a_dense_all = {};
a_cuted = {};
dists_cuted = {};
dists_all = {};
indxs_cuted = {};
integral_ratios_vec = [];
integrals_all = {};
W = 5;
MM = 4;
for i = qs
    x0 = ptfs(:,iter);
    [as, as_dense, dists_dense, samples] = compute_a_sequence(sigma, mu, i, max_vol, M, x0, N, R, r);
    %as_all(iter, :) = as;
    a_dense_all{iter} = as;
    dists_all{iter} = dists_dense;
    %integral_ratios{iter} = compute_integral_ratios_fixed_q_all(sigma, mu, i, as, x0, N, R, r);
    int_ratios = compute_integral_ratios_fixed_q_all_opt(sigma, mu, i, as, x0, N, R, r, samples);
    integrals_all{iter} = int_ratios;
    [as, dists, indxs] = cut_as_3(as, dists_dense, MM);
    a_cuted{iter} = as;
    dists_cuted{iter} = dists;
    int_ratios = int_ratios(indxs);
    indxs_cuted{iter} = indxs;
    integral_ratios{iter} = int_ratios;
    integral_ratios_vec = [integral_ratios_vec int_ratios];
    iter = iter + 1
    length(qs)
end
integral_ratios_vec(isnan(integral_ratios_vec))=0;
integral_ratios_vec(integral_ratios_vec>1) = 1

risk_fun = @(x) medium_risk(x);
dispersion_fun = @(x) low_dispersion(x);

w = compute_weights(vols, dists_all, indxs_cuted, risk_fun, dispersion_fun);

[Ws, Ts] = compute_multiple_weights(w', N*3);

