function [a_vals, samples] = compute_a_dense_sequence(sigma, mu, q, a_max, x0, N)

    curr_fn=3;
    %curr_its=1;
    
    a_stop = 0;
    it = 1;
    a_vals(it) = a_max;
    
    samples = {};
            
    while (curr_fn>1.5 & a_vals(it) >= a_stop)
        %curr_its = N;
        a = a_vals(it);
        it = it+1;
        
        [new_a, pts] = compute_a_dense_next(sigma, mu, q, a, x0, N);
        samples{it - 1} = pts;
        a_vals(it) = new_a
       
        q_mu = q*(mu'*pts);
    
        Sx = sigma*pts;
        log_xSx = zeros(1,N);
    
        for i=1:N
            log_xSx(i) = pts(:,i)' * Sx(:,i);
        end
        curr_fn = mean(exp((- new_a + a) * (log_xSx - q_mu)))
        if (curr_fn < 1)
            curr_fn = 1 / curr_fn
        end
        %mean(exp(( new_a - a) * (log_xSx - q_mu)))
    end
    
    if a_vals(it)>=a_stop
        samples{length(a_vals)} = get_samples(sigma, mu, a_vals(it), q, N, x0);
        %a_vals=a_vals(1:end-1);
    else
        a_vals=a_vals(1:end-1);
    end

end
