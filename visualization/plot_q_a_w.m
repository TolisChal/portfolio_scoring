a=[1 1 0.7 1.2 1 1];
w=[0.25 0.20 0 0 0.3 0.25];
q=[0.3 0.3 1.5 1.5 20 20];

E1 = Evaluate_points3(X,S,mu,1,0.3);
E3 = Evaluate_mixture3(X,S,mu,a,q,w);
scatter3(X(1,:),X(2,:),X(3,:),25,E1,'.')

scatter3(X(1,:),X(2,:),X(3,:),25,E1,'.')
E2 = Evaluate_points3(X,S,mu,1,1.5);
E15 = Evaluate_points3(X,S,mu,1,20);
subplot(1,2,1)
scatter3(X(1,:),X(2,:),X(3,:),25,E1,'.')
X2 = X+0.5;
X3= X2+0.5;
X4 = X3+0.5;
hold on

scatter3(X2(1,:),X2(2,:),X2(3,:),25,E15,'.')
scatter3(X3(1,:),X3(2,:),X3(3,:),25,E2,'.')


[~,r]=max(E1);
scatter3(X(1,r),X(2,r),X(3,r),'*k')

[~,r]=max(E15);
scatter3(X2(1,r),X2(2,r),X2(3,r),'*k')

[~,r]=max(E2);
scatter3(X3(1,r),X3(2,r),X3(3,r),'*k')

subplot(1,2,2)
scatter3(X(1,:),X(2,:),X(3,:),25,E3,'.')
hold on
[~,r]=max(E3);
scatter3(X(1,r),X(2,r),X(3,r),'*k')
