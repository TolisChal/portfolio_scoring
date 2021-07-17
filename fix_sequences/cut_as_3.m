function [as, dists_out, index_cut] = cut_as_3(as, dists, M)
    
    N = length(dists);
    W = ceil(N/M);
    index_cut = 1:W:N;
    
    as = as(index_cut);
    dists_out = dists(index_cut);
    
    
end



