Rets = Returns_12_coins(461:520, :);
[sigma2, ~] = covCor(Rets);
R2 = mean(Rets);

Rs = mvnrnd(R2,sigma2,10000);

X = sample_from_mixture2(ones(70,1)/70, samples_all, 70000);
scores1 = get_scores(OptMVPtf, Rs, X);
y1=ones(1,length(0:0.001:1))*0;

X = sample_from_mixture2(Ws(10,:), samples_all, 70000);
scores2 = get_scores(OptMVPtf, Rs, X);
y2=ones(1,length(0:0.001:1))*0.1;

X = sample_from_mixture2(Ws(20,:), samples_all, 70000);
scores3 = get_scores(OptMVPtf, Rs, X);
y3=ones(1,length(0:0.001:1))*0.2;

X = sample_from_mixture2(Ws(30,:), samples_all, 70000);
scores4 = get_scores(OptMVPtf, Rs, X);
y4=ones(1,length(0:0.001:1))*0.3;

X = sample_from_mixture2(Ws(40,:), samples_all, 70000);
scores5 = get_scores(OptMVPtf, Rs, X);
y5=ones(1,length(0:0.001:1))*0.4;

X = sample_from_mixture2(Ws(50,:), samples_all, 70000);
scores6 = get_scores(OptMVPtf, Rs, X);
y6=ones(1,length(0:0.001:1))*0.5;

X = sample_from_mixture2(Ws(60,:), samples_all, 70000);
scores7 = get_scores(OptMVPtf, Rs, X);
y7=ones(1,length(0:0.001:1))*0.6;

X = sample_from_mixture2(Ws(70,:), samples_all, 70000);
scores8 = get_scores(OptMVPtf, Rs, X);
y8=ones(1,length(0:0.001:1))*0.7;

X = sample_from_mixture2(Ws(80,:), samples_all, 70000);
scores9 = get_scores(OptMVPtf, Rs, X);
y9=ones(1,length(0:0.001:1))*0.8;

X = sample_from_mixture2(Ws(90,:), samples_all, 70000);
scores10 = get_scores(OptMVPtf, Rs, X);
y10=ones(1,length(0:0.001:1))*0.8;

%hist(scores,100)

scores1(scores1==0)=0.0000001;
scores1(scores1==1)=0.9999999;

scores2(scores2==0)=0.0000001;
scores2(scores2==1)=0.9999999;

scores3(scores3==0)=0.0000001;
scores3(scores3==1)=0.9999999;

scores4(scores4==0)=0.0000001;
scores4(scores4==1)=0.9999999;

scores5(scores5==0)=0.0000001;
scores5(scores5==1)=0.9999999;

scores6(scores6==0)=0.0000001;
scores6(scores6==1)=0.9999999;

scores7(scores7==0)=0.0000001;
scores7(scores7==1)=0.9999999;

scores8(scores8==0)=0.0000001;
scores8(scores8==1)=0.9999999;

scores9(scores9==0)=0.0000001;
scores9(scores9==1)=0.9999999;

[PDF1,xi] = ksdensity(scores1, 0:0.001:1, 'Support',[0 1], 'Function', 'pdf');
[PDF2,~] = ksdensity(scores2, 0:0.001:1, 'Support',[0 1], 'Function', 'pdf');
[PDF3,~] = ksdensity(scores3, 0:0.001:1, 'Support',[0 1], 'Function', 'pdf');
[PDF4,~] = ksdensity(scores4, 0:0.001:1, 'Support',[0 1], 'Function', 'pdf');
[PDF5,~] = ksdensity(scores5, 0:0.001:1, 'Support',[0 1], 'Function', 'pdf');
[PDF6,~] = ksdensity(scores6, 0:0.001:1, 'Support',[0 1], 'Function', 'pdf');
[PDF7,~] = ksdensity(scores7, 0:0.001:1, 'Support',[0 1], 'Function', 'pdf');
[PDF8,~] = ksdensity(scores8, 0:0.001:1, 'Support',[0 1], 'Function', 'pdf');
[PDF9,~] = ksdensity(scores9, 0:0.001:1, 'Support',[0 1], 'Function', 'pdf');

figure;
plot3(xi, y1, PDF1);
hold on
plot3(xi, y2, PDF2);
plot3(xi, y3, PDF3);
plot3(xi, y4, PDF4);
plot3(xi, y5, PDF5);
plot3(xi, y6, PDF6);
plot3(xi, y7, PDF7);
plot3(xi, y8, PDF8);
plot3(xi, y9, PDF9);
title(sprintf('Distribution of the score of MV Ptf'));
xlabel('Score');
legend('i = 0', 'i=50', 'i=100', 'i=150', 'i=200', 'i=250', 'i=300', 'i=350', 'i=400');

