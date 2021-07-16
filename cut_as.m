function [as, dists, indxs] = cut_as(as, dists)
    
    M = length(dists);
    
    dists_out = dists - 1;
    dists_out(1) = 0;
    
    index_stop = 1;
    val = dists_out(2);
    for j = 3:M
        if(val / dists_out(j)
    
    i = 1;
    for j = 1:length(indxs)
        dists_out(i) = sum(dists()
    end
    
    as = as(indxs);
    dists = dists(indxs);

end



