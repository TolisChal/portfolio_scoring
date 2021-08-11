function X = sample_from_mixture2(w, samples_all, N)

    n = size(samples_all{1},1);
    
    %M = length(w);
    %num_mix = 0;
    X = zeros(n, N);
    
    for i=1:N
        
        r = rand;
        pos = find(r<cumsum(w));
        pos = pos(1);
        
        s = samples_all{pos};
        M = size(s, 2);
        j = randi([1 M],1,1);
        X(:,i) = s(:,j);
    end

end

