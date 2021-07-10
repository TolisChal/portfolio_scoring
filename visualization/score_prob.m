y0 = Sampling_simplex(3,1,'RM');
mu=randn(1,3);
S = randn(7,3);
S=S'*S;
N = 500000;

X = Sampling_simplex(3,N,'RM');

X2 = X - repmat(y0,[1 N]);

m2 = mu*X2;
Sx = S*X2;
s2 = zeros(1,N);
evals = zeros(1,N);

for i=1:N
    s2(i) = X2(:,i)' * Sx(:,i);
    evals(i) = 1 - 0.5*(1 + erf((0-m2(i))/s2(i)));
end

scatter3(X(1,:),X(2,:),X(3,:),25,evals,'.')
hold on
scatter3(y0(1), y0(2), y0(3), '*k')
[~,r]=max(evals);
scatter3(X(1,r),X(2,r),X(3,r),'*r')
c = colorbar
sum(evals==1)

figure(2)

x=0:0.01:1;
for i=1:length(x)
    y(i) = sum(evals<=x(i))/N;
end
for i=1:(length(y)-1)
   dy(i) = (y(i+1)-y(i))/0.01; 
end
plot(x,y,'.')


figure(3)
plot(x(1:(length(x)-1)),dy,'.')

figure(4)
plot(y,x,'.')

figure(5)

for i=1:(length(y)-1)
   dx(i) = (x(i+1)-x(i))/(y(i+1)-y(i)); 
end
plot(y(1:(length(y)-1)),dx,'.')