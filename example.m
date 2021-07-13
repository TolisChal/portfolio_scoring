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
[qs, ptfs] = compute_q_sequence(sigma, mu, min_vol, max_vol, 5);

q=qs(2);
a=2;

X = Sampling_simplex(3, 10000, 'RM');

vals = Evaluate_density(X, sigma, mu, a, q);

scatter3(X(1,:),X(2,:),X(3,:),25,vals,'.')

x0 = ptfs(:,2);

[as, as_dense, samples] = compute_a_sequence(sigma, mu, q, max_vol, M, x0, N);

Ptf = Sampling_simplex(3,1, 'RM');

R = randn(1,3);
r = R*Ptf;

integral_ratios = compute_integral_ratios_fixed_q(sigma, mu, q, as, as_dense, x0, N, R, r);

for i = as
    figure
    Y = get_samples(sigma, mu, i, q, 5000, x0);

    subplot(1,2,1)
    plot3(Y(1,:), Y(2,:), Y(3,:), '.')
    xlim([0 1])
    ylim([0 1])
    zlim([0 1])
    
    Y = get_samples_restricted(sigma, mu, a, q, N, R, r);

    subplot(1,2,2)
    plot3(Y(1,:), Y(2,:), Y(3,:), '.')
    xlim([0 1])
    ylim([0 1])
    zlim([0 1])
end




