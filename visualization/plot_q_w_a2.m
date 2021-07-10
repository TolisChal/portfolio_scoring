scatter3(X(1,:),X(2,:),X(3,:),25,E1,'.')
E1 = Evaluate_points3(X,S,mu,1,0.3);
scatter3(X(1,:),X(2,:),X(3,:),25,E1,'.')
E2 = Evaluate_points3(X,S,mu,1,2);
E15 = Evaluate_points3(X,S,mu,1,1);
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
c=colorbar
subplot(1,2,2)
Xit = X;
for i=1:3
    i
a=rand(1,6)*2
w=Sampling_simplex(6,1,'RM')'
q=rand*2;
q=[q q];
q2=rand*2;
q=[q q2 q2];
q2=rand*2;
q=[q q2 q2]

E3 = Evaluate_mixture2(X,S,mu,a,q,w);


scatter3(Xit(1,:),Xit(2,:),Xit(3,:),25,E3,'.')

hold on
[~,r]=max(E3);
scatter3(Xit(1,r),Xit(2,r),Xit(3,r),'*k')
Xit = Xit+0.5;
end
c=colorbar
