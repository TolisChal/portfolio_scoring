function [as, dists_out, index_cut] = cut_as_2(as, dists, W)
    
    M = length(dists);
    index_cut = 1:W:M;
    
    as = as(index_cut);
    dists_out = dists(index_cut);
    
    return
    
    dists_out = dists - 1;
    dists_out(1) = 0
    
    index_stop = M;
    val = dists_out(2);
    for j = 3:M
        dists_out(j) / val
        if(dists_out(j) / val < 0.1)
            index_stop = j;
            break
        end
    end
    
    as = as(1:index_stop);
end



