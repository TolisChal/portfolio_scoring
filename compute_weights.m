function w = compute_weights(as, qs, dists_dense, risk_fun, dispersion_fun, q_min, q_max)
    
    qs = [q_min qs q_max];
    M1 = length(qs);
    q = 0:(1/M1):1;
    q_vals = risk_fun(q);
    
    q_vals = q_vals(2:(M1-1));
    M1 = M1 -2;
    
    w= [];
    
    for i = 1:M1
        a_seq = as{i};
        dists = dists_dense{i};
        interval_length = sum(dists);
        cum_sum_dists = cumsum(dists);
        
        a_vals_int = cum_sum_dists / interval_length;
        a_vals_int = sort(1 - a_vals_int);
        
        a_vals = dispersion_fun(a_vals_int);
        
        wi = q_vals(i);
        w = [w wi*a_vals];
    end
    
    w = w / norm(w);
end