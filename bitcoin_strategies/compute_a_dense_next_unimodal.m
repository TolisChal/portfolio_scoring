function [new_a, pts] = compute_a_dense_next_unimodal(sigma, y0, a, N)

    n = length(y0);
    ratio = 1 - 1/n;
    C = 4;
    last_a = a;
            
    done = 0;
    k = 1;
    pts = get_samples_unimodal(sigma, y0, last_a, N, y0);% (sigma, mu, last_a, q, N, x0);
    
    X = pts - repmat(y0, [1 size(pts, 2)]);
    Sx = sigma*X;
    log_xSx = zeros(1,size(X, 2));
    
    for i=1:size(X, 2)
        log_xSx(i) = X(:,i)' * Sx(:,i);
    end
    
    last_ratio = 0.1
    %fn = zeros(its,1);
    while ~done
        a = last_a * ratio^k

        fn = exp((- a + last_a) * (log_xSx));
               
        
        (var(fn)/mean(fn))/mean(fn)
        var(fn)/(mean(fn)^2)
        %mean(fn)/last_ratio
        %mean(fn)
        ratio_iter = mean(fn)/last_ratio
        if (ratio_iter < 1)
            ratio_iter = 1/ratio_iter
        end
        done=1;
        if ((var(fn)/mean(fn))/mean(fn)<C || ratio_iter<1+1e-5)
            if k~=1
                
                k = k/2;
            end
            done = 1;
        else
            k = 2*k;
        end
        last_ratio = mean(fn);
    end
            
    new_a = last_a * ratio^k;
end

