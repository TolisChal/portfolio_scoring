function Copula = compute_copula(sigma, R, m, X)
    
    N = size(X,2);
    
    RndPtfRet_P1 = R * X;
    Ex = sigma * X;
    
    Vol = zeros(N, 1);
    for j=1:N
        Vol(j) = X(:, j)' * Ex(:, j);
    end

    [~, I1] = sort(RndPtfRet_P1,'ascend');
    [~, I2] = sort(Vol,'ascend');

    for j=1:m
        range = ((j-1)*floor(N/m)+1:j*floor(N/m));
        ScoreP1(I1(range))=j;
        Volatility(I2(range))=j;
    end

    Copula = zeros(m, m);
    for k=1:m
        for j=1:m
            Copula(k,j) = sum((ScoreP1==k).*(Volatility==j));
        end
    end
    %sum(sum(Copula))
    Copula=Copula/N;
    %copulas{i} = Copula;

end


