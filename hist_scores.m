Rets = Returns_12_coins(461:520, :);
[sigma2, ~] = covCor(Rets);
R2 = mean(Rets);

Rs = mvnrnd(R2,sigma2,10000);

X = sample_from_mixture2(Ws(1,:), samples_all, 70000);

scores = get_scores(OptMVPtf, Rs, X);

%hist(scores,100)

scores(scores==0)=0.0000001;
scores(scores==1)=0.9999999;
[PDF1,xi] = ksdensity(scores, 0:0.001:1, 'Support',[0 1], 'Function', 'pdf');

plot(xi, PDF1);
hold on
title(sprintf('Distribution of the score'));
xlabel('Score');
legend('Ptf Avg. volatility');

