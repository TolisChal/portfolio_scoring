addpath('./fix_sequences')
addpath('./mcmc_diagnostics')
addpath('./compute_integrals')

d = 3;

N = 10000;

mu = rand(d,1);
sigma = rand(d);
sigma = sigma'*sigma + eye(d)*3;
sigma = sigma / (0.34*max(max(sigma)));

[min_vol, max_vol] = compute_min_max_volatility(sigma);
[qs, ptfs] = compute_q_sequence(sigma, mu, min_vol, max_vol, 5);



q=qs(1);
a=2;

X = Sampling_simplex(3, 10000, 'RM');

vals = Evaluate_density(X, sigma, mu, a, q);

scatter3(X(1,:),X(2,:),X(3,:),25,vals,'.')

x0 = ptfs(:,1);

[as, samples] = compute_a_sequence(sigma, mu, q, M, x0, N);

for i = as
    Y = get_samples(sigma, mu, i, q, 5000, x0);

    figure
    plot3(Y(1,:), Y(2,:), Y(3,:), '.')
    xlim([0 1])
    ylim([0 1])
    zlim([0 1])
end
%N = 5000;
%as = compute_a_sequence(sigma, mu, q, 4, max_vol, x0, N);


