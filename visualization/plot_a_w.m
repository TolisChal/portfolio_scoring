X = Sampling_simplex(3,150000,'RM');
X2 = X+0.5;
X3 = X2+0.5;
X4 = X3+0.5;
a=[0.5 1 1.5];
w=[0.27 0.4 0.33];
E1 = Evaluate_points(X,S,mu,1);
E2 = Evaluate_points(X,S,mu,1.5);
E3 = Evaluate_points(X,S,mu,2);
E4 = Evaluate_mixture1(X,S,mu,a,w);
subplot(1,2,1)
scatter3(X(1,:),X(2,:),X(3,:),25,E1,'.')
hold on
scatter3(X2(1,:),X2(2,:),X2(3,:),25,E2,'.')
scatter3(X3(1,:),X3(2,:),X3(3,:),25,E3,'.')


[~,r]=max(E1);
scatter3(X(1,r),X(2,r),X(3,r),'*k')

[~,r]=max(E2);
scatter3(X2(1,r),X2(2,r),X2(3,r),'*k')

[~,r]=max(E3);
scatter3(X3(1,r),X3(2,r),X3(3,r),'*k')
c=colorbar

subplot(1,2,2)
scatter3(X(1,:),X(2,:),X(3,:),25,E4,'.')
hold on
[~,r]=max(E4);
scatter3(X(1,r),X(2,r),X(3,r),'*k')
c=colorbar
%----------------------------------------------------------%

