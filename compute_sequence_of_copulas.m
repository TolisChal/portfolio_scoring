Rets = Returns_12_coins(461:520, :);
[sigma2, ~] = covCor(Rets);
R2 = (prod(1 + Rets) - 1);
copulas_low = {};

X = sample_from_mixture2(ones(70,1)/70, samples_all, 70000);
Copula = compute_copula(sigma2, R2, 50, X);
copulas_low{1} = Copula;

WW = 50;
iter = 2;
for i=50:WW:400
    X = sample_from_mixture2(Ws3(i,:), samples_all, 70000);
    Copula = compute_copula(sigma2, R2, 50, X);
    copulas_low{iter} = Copula;
    iter = iter + 1;
end

X = sample_from_mixture2(Ws3(422,:), samples_all, 70000);
Copula = compute_copula(sigma2, R2, 50, X);
copulas_low{iter} = Copula;

indis_low = comp_indicators_mask(copulas_low);
for j=2:(length(copulas_low)-1)
    figure
    surf(copulas_low{j})
    title(strcat('Indicator = ',num2str(indis_low(j))))
    xlabel('Portfolio return')
    ylabel('Portfolio volatility')
    set(gca,'ZTick',[])
    xticks([0 50])
    xticklabels({'0', '1'})
    yticks([0 50])
    yticklabels({'0', '1'})
    grid on
end


