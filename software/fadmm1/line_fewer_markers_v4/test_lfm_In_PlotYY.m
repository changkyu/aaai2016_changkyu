clear all; close all;

figure(10); hold on; set(gca,'FontSize',16); set(gca,'FontName','Times'); set(gcf,'Color',[1,1,1]);
x=(1:10);
y1=x.^2;y2=x.^3;

f1 = @(x,y)line_fewer_markers(x,y,5,'ro');
f2 = @(x,y)line_fewer_markers(x,y,8,'-.bs','Spacing','curve','markerfacecolor','g');%,'LegendLine','off');
[yH,lh1,lh2] = plotyy(x,y1,x,y2,f1,f2)

xlabel('x');
ylabel(yH(1),'x^2')
ylabel(yH(2),'x^3')
legend([lh1,lh2],{'curve 1','curve 2'},'location','northwest')
set(gcf,'position',[100 100 900 600]);