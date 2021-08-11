nassets = 3;

df = 80;
sigma =  wishrnd(eye(nassets),df)/df;
mu = randn(nassets,1);
asset_returns = randn(1,nassets);

dispersion_fun = @(x) low_dispersion(x);
risk_fun = @(x) medium_risk(x);

[parametric_score, Ws, bias_vector, mean_score, a_sequence, q_sequence, volatility_sequence] = compute_scores(sigma, mu, asset_returns, 3, 4, dispersion_fun, dispersion_fun);
