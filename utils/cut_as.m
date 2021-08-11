function [as, dists_out, index_stop] = cut_as(as, dists)
    
    M = length(dists);
    
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



