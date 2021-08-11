addpath('./fix_sequences')
addpath('./mcmc_diagnostics')
addpath('./compute_integrals')

close all

d = 3;
M = 5;

N = 10000;

mu = rand(d,1);
sigma = rand(d);
sigma = sigma'*sigma + eye(d)*3;
sigma = sigma / (0.34*max(max(sigma)));

[min_vol, max_vol] = compute_min_max_volatility(sigma);
[qs, ptfs, q_min, q_max, vols] = compute_q_sequence(sigma, mu, min_vol, max_vol, 5);

q=qs(2);
a=2;

X = Sampling_simplex(3, 10000, 'RM');

vals = Evaluate_density(X, sigma, mu, a, q);

scatter3(X(1,:),X(2,:),X(3,:),25,vals,'.')


Ptf = Sampling_simplex(3,1, 'RM');

R = randn(1,3);
r = R * Ptf;

integral_ratios = {};
iter = 1;
%as_all = zeros(5, 5);
a_dense_all = {};
dists_all = {};
integral_ratios_vec = [];
for i = qs
    x0 = ptfs(:,iter);
    [as, as_dense, dists_dense, samples] = compute_a_sequence(sigma, mu, i, max_vol, M, x0, N);
    %as_all(iter, :) = as;
    a_dense_all{iter} = as;
    dists_all{iter} = dists_dense;
    %integral_ratios{iter} = compute_integral_ratios_fixed_q_all(sigma, mu, i, as, x0, N, R, r);
    integral_ratios{iter} = compute_integral_ratios_fixed_q_all_opt(sigma, mu, i, as, x0, N, R, r, samples);
    integral_ratios_vec = [integral_ratios_vec integral_ratios{iter}];
    iter = iter + 1
    length(qs)
end
integral_ratios_vec(isnan(integral_ratios_vec))=0

risk_fun = @(x) medium_risk(x);
dispersion_fun = @(x) low_dispersion(x);

w = compute_weights(a_dense_all, qs, dists_all, risk_fun, dispersion_fun, q_min, q_max);

X = sample_exponential(w', 50, N);
x_c = mean(X,2);
%for i = as
%    figure
%    Y = get_samples(sigma, mu, i, q, 5000, x0);

%    subplot(1,2,1)
%    plot3(Y(1,:), Y(2,:), Y(3,:), '.')
%    xlim([0 1])
%    ylim([0 1])
%    zlim([0 1])
    
%    Y = get_samples_restricted(sigma, mu, a, q, N, R, r);

%    subplot(1,2,2)
%    plot3(Y(1,:), Y(2,:), Y(3,:), '.')
%    xlim([0 1])
%    ylim([0 1])
%    zlim([0 1])
%end




