x1=0:0.001:1;
y1 = exp(x1)-0.8;

x2 = 0:0.001:1;
y2 = 8*(exp(-(x2 - 0.5).^2)-exp(-0.25)) + 0.2;

plot(x1, y1, 'r.')
hold on
plot(x2,y2, 'k.')