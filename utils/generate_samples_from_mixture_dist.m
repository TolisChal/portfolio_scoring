samples_all = {};
N = 20000;
M1 = 3;
M2 = 4;

for i=1:length(qs)
    as = a_cuted{i};
    x0 = ptfs(:,i);
    q = qs(i);
    samples = {};
    for j=1:length(as)
        disp([i j])
        a=as(j);
        X = get_samples(sigma, mu, a, q, N, x0);
        samples{j} = X;
    end
    samples_all{i} = samples;
end

save('simple_example_12_strategies_samples_high_risk.mat','samples_all')
