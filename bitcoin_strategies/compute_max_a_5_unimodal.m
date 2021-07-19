function a_max = compute_max_a_5_unimodal(sigma, y0, N, a0_max, R, r)
    
    n = length(y0);
    ratio = 1 + (1/sqrt(n));
    
    
    a_max = a0_max;
    

    while(true)

        a_max = a_max*ratio
        
        X = get_samples_unimodal(sigma, y0, a_max, N, y0); %(sigma, mu, a_max, q, N, x0);
        
        vecc = sqrt(sum((X - repmat(y0,[1, size(X,2)])).^2,1));
        length(vecc)
        rad = max(vecc)
        
        sc=R*X;
        sc = sum(sc<r) / size(X,2)
        %left
        if (sc>0.95 & rad < 0.11)
            return
        elseif (sc<0.05 & rad < 0.11)
            return
        end
        
    end
    
end