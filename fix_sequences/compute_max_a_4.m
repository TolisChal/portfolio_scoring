function a_max = compute_max_a_4(sigma, mu, q, N, a0_max, x0, R, r)
    
    left = false;
    
    if(R*x0 < r)
        left = true;
    end
    
    n = length(mu);
    ratio = 1 + (1/sqrt(n));
    
    
    a_max = a0_max;
    

    while(true)

        a_max = a_max*ratio
        
        X = get_samples(sigma, mu, a_max, q, N, x0);
        
        vecc = sqrt(sum((X - repmat(x0,[1, size(X,2)])).^2,1));
        length(vecc)
        max(vecc)
        
        sc=R*X;
        sc = sum(sc<r) / size(X,2)
        left
        if (left & sc>0.95)
            return
        elseif (~left & sc<0.05)
            return
        end
        
    end
    
end