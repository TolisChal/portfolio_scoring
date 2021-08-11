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

[as, samples] = compute_a_sequence(sigma, mu, q, max_vol, M, x0, N);

Ptf = Sampling_simplex(3,1, 'RM');

a = as(4);

R = randn(1,3);
r = R*Ptf;

Y = get_samples(sigma, mu, a, q, N, x0);

figure
plot3(Y(1,:), Y(2,:), Y(3,:), '.')
xlim([0 1])
ylim([0 1])
zlim([0 1])

Y2 = get_samples_restricted(sigma, mu, a, q, N, R, r);

figure
plot3(Y2(1,:), Y2(2,:), Y2(3,:), '.')
xlim([0 1])
ylim([0 1])
zlim([0 1])

Y3 = billiard_walk_low_dim(R, r, 2, N);

figure
plot3(Y3(1,:), Y3(2,:), Y3(3,:), '.')
xlim([0 1])
ylim([0 1])
zlim([0 1])

X = Sampling_simplex(3, 10000, 'RM');

sc = R * X;
sc = sum(sc < r) / N

int_ratio1 = estimate_integral_ratio(sigma, mu, a, q, N, x0, R, r);
int_ratio2 = estimate_integral_ratio_naive(sigma, mu, a, q, N, R, r)

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


