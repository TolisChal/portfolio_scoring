
btc_ptf = zeros(12,1);
btc_ptf(1) = 1;
r = R*ptf_erc;
erc_scores_high=zeros(size(Ws,1),1);
for i=1:size(Ws,1)
    i
    X = sample_from_mixture2(Ws(i,:), samples_all, 70000);
    sc = R*X;
    erc_scores_high(i) = sum(sc<r) / length(sc);
end

%plot(1:size(Ws,1),scores,'r.')
%xlabal('i-th iteration ()')

ll = min([length(OptMV_scores_high) length(OptMV_scores_medium) length(OptMV_scores_low)]);


p=plot(1:ll, OptMV_scores_high(1:ll), '-r.', 1:ll, OptMV_scores_medium(1:ll), '-m.', 1:ll, OptMV_scores_low(1:ll), '-b.');
p(1).MarkerSize = 14;
p(2).MarkerSize = 14;
p(3).MarkerSize = 14;
h=legend('Plot D (High risk)', 'Plot B (medium risk)', 'Plot C (low risk)', 'Location','northwest');
h.FontSize = 20;
grid on
xlabel('iteration','fontsize',20)
ylabel('score','fontsize',20)
ylim([0 1])
title('Parametric score of MV Ptf','fontsize',24)


ll = min([length(erc_scores_high) length(erc_scores_medium) length(erc_scores_low)]);

figure;
p=plot(1:ll, erc_scores_high(1:ll), '-r.', 1:ll, erc_scores_medium(1:ll), '-m.', 1:ll, erc_scores_low(1:ll), '-b.');
p(1).MarkerSize = 14;
p(2).MarkerSize = 14;
p(3).MarkerSize = 14;
h=legend('Plot D (High risk)', 'Plot B (medium risk)', 'Plot C (low risk)', 'Location','northwest');
h.FontSize = 20;
grid on
xlabel('iteration','fontsize',20)
ylabel('score','fontsize',20)
ylim([0 1])
title('Parametric score of ETC Ptf','fontsize',24)



ll = min([length(btc_scores_high) length(btc_scores_medium) length(btc_scores_low)]);

figure;
p=plot(1:ll, btc_scores_high(1:ll), '-r.', 1:ll, btc_scores_medium(1:ll), '-m.', 1:ll, btc_scores_low(1:ll), '-b.');
p(1).MarkerSize = 14;
p(2).MarkerSize = 14;
p(3).MarkerSize = 14;
h=legend('Plot D (High risk)', 'Plot B (medium risk)', 'Plot C (low risk)', 'Location','northeast');
h.FontSize = 20;
grid on
xlabel('iteration','fontsize',20)
ylabel('score','fontsize',20)
ylim([0 1])
title('Parametric score of BTC Ptf','fontsize',24)
