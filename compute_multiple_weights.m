function [Ws, Ts] = compute_multiple_weights(w, N)
    
    n = length(w);
    T_max = sqrt(n)/100;
    
    T_max = compute_a_max_exponential(w, N, T_max);
    %T_min = compute_a_min_exponential(w, N, T_min);
    
    T = T_max
    ratio = 1 - (1/n);
    Ws = [];
    Ts = [];
    
    opt_ptf = zeros(n,1);
    opt_ptf(find(w==max(w))) = 1
    
    while(true)
        
        T = T * ratio;
        X = sample_exponential(-w, T, N); %, A_add, b_add)
        xc = mean(X,2)
        
        Ws = [Ws; xc'];
        Ts = [Ts T];
        
        sum(abs(xc - opt_ptf))
        if (sum(abs(xc - opt_ptf)) < 0.05)
            return
        end
        
    end


end