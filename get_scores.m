function scores = get_scores(x0, Rets, X)

    N = size(Rets,1);
    M = size(X,2);
    scores = zeros(1, N);
    
    for i=1:N
        i
        R = Rets(i,:);
        r = R*x0;
        sc = R*X;
        
        scores(i) = sum(sc<r)/M;
    end

end