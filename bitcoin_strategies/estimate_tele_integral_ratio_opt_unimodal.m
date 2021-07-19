function [int_ratio, const_int_num, const_int_den] = estimate_tele_integral_ratio_opt_unimodal(sigma, y0, as_merged, current_pos, to_pos, R, r, const_int_num, const_int_den, volume_ratio, N, samples)

    n = length(y0)
    k = length(as_merged);
    
    for i = current_pos:to_pos

        %X1 = get_samples(sigma, mu, as_merged(i), q, N, x0);
        X1 = samples{k + 1 - i};
        X2 = get_samples_restricted_unimodal(sigma, y0, as_merged(i), N, R, r);
        
        X1 = X1 - repmat(y0, [1, size(X1, 2)]);
        Sx = sigma*X1;
        log_xSx = zeros(1,size(X1, 2));
        for j=1:size(X1,2)
            log_xSx(j) = X1(:,j)' * Sx(:,j);
        end
        const_int_den = const_int_den * mean(exp((-as_merged(i+1) + as_merged(i))*(log_xSx)));
        
        X2 = X2 - repmat(y0, [1, size(X2, 2)]);
        Sx = sigma*X2;
        log_xSx = zeros(1,size(X2,2));
        for j=1:size(X2,2)
            log_xSx(j) = X2(:,j)' * Sx(:,j);
        end
        const_int_num = const_int_num * mean(exp((-as_merged(i+1) + as_merged(i))*(log_xSx)));
    end
    
    int_ratio = volume_ratio * (const_int_num / const_int_den);

end

