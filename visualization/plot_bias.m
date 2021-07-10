c = randn(1,3);
c=c/norm(c)

E = exp(0.01*c*X);

scatter3(X(1,:),X(2,:),X(3,:),25,E,'.')
hold on
scatter3(0,0,1,'*k')
scatter3(1,0,0,'*r')
scatter3(0,1,0,'*g')

xlabel('x label')
ylabel('y label')
zlabel('z label')

