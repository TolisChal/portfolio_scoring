function [parametric_score, Ws, bias_vector, mean_score, a_sequence, q_sequence, volatility_sequence] = compute_scores(sigma, mu, asset_returns, Ptf, num_risk_levels, num_dispersion_levels, risk_behavioral_function, dispersion_behavioral_function)

    N = 10000;
    nassets = length(mu);

    [min_vol, max_vol] = compute_min_max_volatility(sigma);
    [qs, ptfs, q_min, q_max, vols] = compute_q_sequence(sigma, mu, min_vol, max_vol, num_risk_levels);

    %[OptMVPtf, avg_vol] = compute_2avg_volatility_ptf(sigma, mu);
    OptMVPtf = Ptf;
    %OptMVPtf
    R = asset_returns;
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
    samples_all = {};
    W = 5;
    MM = num_dispersion_levels;
    
    for i = qs
        x0 = ptfs(:,iter);
        [as, as_dense, dists_dense, samples] = compute_a_sequence(sigma, mu, i, max_vol, MM, x0, N, R, r);
        %as
        %as_all(iter, :) = as;
        a_dense_all{iter} = as;
        dists_all{iter} = dists_dense;
        %integral_ratios{iter} = compute_integral_ratios_fixed_q_all(sigma, mu, i, as, x0, N, R, r);
        int_ratios = compute_integral_ratios_fixed_q_all_opt(sigma, mu, i, as, x0, N, R, r, samples);
        integrals_all{iter} = int_ratios;
        [as, dists, indxs] = cut_as_3(as, dists_dense, MM);
        for j=1:length(indxs)
            samples_all{length(samples_all)+1} = samples{indxs(j)};
        end
        a_cuted{iter} = as;
        dists_cuted{iter} = dists;
        int_ratios = int_ratios(indxs);
        indxs_cuted{iter} = indxs;
        integral_ratios{iter} = int_ratios;
        integral_ratios_vec = [integral_ratios_vec int_ratios];
        iter = iter + 1;
        %length(qs)
    end
    %integral_ratios_vec(isnan(integral_ratios_vec))=0;
    integral_ratios_vec(integral_ratios_vec>1) = 1;

    %risk_fun = @(x) medium_risk(x);
    %dispersion_fun = @(x) low_dispersion(x);
    %low_risk_fun = @(x) low_risk(x);

    w = compute_weights(vols, dists_all, indxs_cuted, risk_behavioral_function, dispersion_behavioral_function);

    [Ws, ~] = compute_multiple_weights(w', N*2);

    parametric_score = Ws * integral_ratios_vec';
    bias_vector = w';
    mean_score = (1/length(w))*ones(1,length(w)) * integral_ratios_vec';
    a_sequence = a_cuted;
    q_sequence = qs;
    volatility_sequence = vols(2:(end-1));
end


