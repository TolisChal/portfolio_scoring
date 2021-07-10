d=10;
y0 = Sampling_simplex(d,1,'RM');
mu=randn(1,d);
S = randn(2*d,d);
S=S'*S;
N = 500000;

evals12=[];

X = Sampling_simplex(d,N,'RM');

for i =1:4
    i
X2 = X - repmat(y0,[1 N]);

m2 = mu*X2;
Sx = S*X2;
s2 = zeros(1,N);
evals = zeros(1,N);

for j=1:N
    s2(j) = X2(:,j)' * Sx(:,j);
    evals(j) = 1 - 0.5*(1 + erf((0-m2(j))/s2(j)));
end

%scatter3(X(1,:),X(2,:),X(3,:),25,evals,'.')
%hold on
%scatter3(y0(1), y0(2), y0(3), '*k')
[~,r]=max(evals);
%scatter3(X(1,r),X(2,r),X(3,r),'*r')
%c = colorbar
%sum(evals==1)

%figure

x=0:0.01:1;
for j=1:length(x)
    y(j) = sum(evals<=x(j))/N;
end
for j=j:(length(y)-1)
   dy(j) = (y(j+1)-y(j))/0.01; 
end
%plot(x,y,'.')


%figure
%plot(x(1:(length(x)-1)),dy,'.')

pos = sum(evals<=0.5)/N
RR = mu*y0
VAR  = sqrt(y0'*S*y0)

ratio = RR/VAR

XX = mvnrnd(mu,S,1000000);
for j=1:1000000
    evals2(j) = Ali73(XX(j,:), XX(j,:)*y0);
end
if(i==1 || i==4)
    jj=1;
    if (i==4)
        jj=2;
    end
    evals12 = [evals12; evals2];
subplot(1,2,jj)
%figure
hist(evals2,100);
hists(i,:) = hist(evals2,100);
ylim([0 50000])

%[PDF2,x2] = ksdensity(evals2, min(evals2):0.001:max(evals2), 'Function', 'pdf');
end

y0 = X(:,r);
end

for i=1:2
    subplot(1,2,i)
    ylim([0 max(max(hists))])
end

[PDF1,x1] = ksdensity(evals12(1,:), min(evals12(1,:)):0.001:max(evals12(1,:)), 'Function', 'pdf');
figure;
plot(x1, PDF1, 'LineWidth', 4);
title(sprintf('Distribution of the score'),'FontSize', 26);
xlabel('Score','FontSize', 30);

[PDF2,x2] = ksdensity(evals12(2,:), min(evals12(2,:)):0.001:max(evals12(2,:)), 'Function', 'pdf');
figure;
plot(x2, PDF2, 'LineWidth', 4);
title(sprintf('Distribution of the score'),'FontSize', 26);
xlabel('Score','FontSize', 30);
%figure(4)
%plot(y,x,'.')

%figure(5)

%for i=1:(length(y)-1)
%    dx(i) = (x(i+1)-x(i))/(y(i+1)-y(i));
%end
%plot(y(1:(length(y)-1)),dx,'.')