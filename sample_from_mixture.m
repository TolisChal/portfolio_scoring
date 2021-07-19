function X = sample_from_mixture(w, samples, N)

    n = length(w);
    
    M = length(samples);
    num_mix = 0;
    samples_all = {};
    
    counter = 1;
    for i = 1:M
        s = samples{i};
        MM = length(s);
        for j = 1:MM
            samples_all{counter} = s{j};
            counter = counter + 1;
        end
    end
    %length(samples_all)
    samples = [];
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

