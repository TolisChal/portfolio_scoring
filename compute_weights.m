function w = compute_weights(vols, dists_dense, indx_cuted, risk_fun, dispersion_fun)
    
    vols_vals = (vols - vols(1))*(1 ./(vols(end) - vols(1)));
    q_vals = risk_fun(vols_vals);
    
    q_vals = q_vals(2:(end-1));
    
    M1 = length(q_vals);
    
    w= [];
    
    for i = 1:M1
        %a_seq = as{i};
        dists = dists_dense{i};
        dists = dists -1;
        dists(1) = 0;
        indx = indx_cuted{i};
        k = length(dists);
        
        sum_cum_dists = cumsum(dists);
        dists2 = (sum_cum_dists - sum_cum_dists(1))*(1 ./(sum_cum_dists(end) - sum_cum_dists(1)));
        dists2 = dists2(indx);
        a_vals_int = sort(1-dists2)
        
        %a_vals_int = a_vals_int(k+1-indx)
        
        a_vals = dispersion_fun(a_vals_int);
        
        wi = q_vals(i);
        vec = wi*a_vals;
        w = [w vec(length(vec):-1:1)];
    end
    
    w = w / sum(w);
end