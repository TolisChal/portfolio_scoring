Rets = Returns_12_coins(461:520, :);
[sigma2, ~] = covCor(Rets);
R2 = mean(Rets);

Rs = mvnrnd(R2,sigma2,10000);

WW= 10;

X = sample_from_mixture2(ones(70,1)/70, samples_all, 70000);
scores1 = get_scores(OptMVPtf, Rs, X);
y1=ones(1,length(0:0.001:1))*0;

scores1(scores1==0)=0.0000001;
scores1(scores1==1)=0.9999999;

[PDF1,xi] = ksdensity(scores1, 0:0.001:1, 'Support',[0 1], 'Function', 'pdf');

ys = y1;
scores_dists = scores1;
pdfs = PDF1;

figure;
plot3(xi, y1, PDF1);
hold on

for i=10:WW:410
    
    i
    X = sample_from_mixture2(Ws3(i,:), samples_all, 70000);
    scores2 = get_scores(OptMVPtf, Rs, X);
    y2=ones(1,length(0:0.001:1))*(i/100);
    
    scores2(scores2==0)=0.0000001;
    scores2(scores2==1)=0.9999999;
    
    [PDF2,~] = ksdensity(scores2, 0:0.001:1, 'Support',[0 1], 'Function', 'pdf');
    
    plot3(xi, y2, PDF2);
    
    ys = [ys; y2];
    scores_dists = [scores_dists; scores2];
    pdfs = [pdfs; PDF2];
    
end

title('Distributions of the score of MV Ptf (Plot B - Medium risk) ', 'fontsize',20);
set(gca,'ZTick',[], 'YTick', [])
xlabel('Score','fontsize',20)
%legend('i = 0', 'i=50', 'i=100', 'i=150', 'i=200', 'i=250', 'i=300', 'i=350', 'i=400');


