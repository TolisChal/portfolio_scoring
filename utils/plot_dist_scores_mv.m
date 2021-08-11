figure
for i=1:2:42
plot3(xi,ys2(i,:),pdfs2(i,:))
hold on
end
title('Distributions of the score of MV Ptf (Plot D - High risk) ', 'fontsize',20);
set(gca,'ZTick',[], 'YTick', [])
xlabel('Score','fontsize',20)
grid on

figure
for i=1:2:32
plot3(xi,ys3(i,:),pdfs3(i,:))
hold on
end
title('Distributions of the score of MV Ptf (Plot B - Medium risk) ', 'fontsize',20);
set(gca,'ZTick',[], 'YTick', [])
xlabel('Score','fontsize',20)
grid on

figure
for i=1:2:42
plot3(xi,ys(i,:),pdfs(i,:))
hold on
end
title('Distributions of the score of MV Ptf (Plot C - Low risk) ', 'fontsize',20);
set(gca,'ZTick',[], 'YTick', [])
xlabel('Score','fontsize',20)
grid on