f = @(x,y) (1-x).^2 + 100*(y-x.^2).^2;
x = linspace(-1.5,1.5); y = linspace(-1,3);
[xx,yy] = meshgrid(x,y); ff = f(xx,yy);
levels = 10:10:300;
LW = 'linewidth'; FS = 'fontsize'; MS = 'markersize';
figure, contour(x,y,ff,levels,LW,1.2), colorbar
axis([-1.5 1.5 -1 3]), axis square, hold on
xlabel('x','Interpreter','Latex');
ylabel('y','Interpreter','Latex');
title([num2str(D),'Rosenbrock´s function in 2D'],'Interpreter','Latex');
set(gca,'TickLabelInterpreter', 'Latex')
set(gcf,'color','w');