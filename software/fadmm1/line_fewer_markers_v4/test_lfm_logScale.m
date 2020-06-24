clear all; close all;

figure; hold on; grid on;
set(gca,'FontSize',16); set(gca,'FontName','Times'); set(gcf,'Color',[1,1,1]);
xlabel('f');ylabel('|H| [dB]');

f = logspace(log10(100),log10(10e6),100);
pS1= -2*pi*1e4;
pS2= -2*pi*1e5;
pS3= -2*pi*1e6;
G1 = -pS1./(pS1+2*pi*j*f);
G2 = -pS2./(pS2+2*pi*j*f);
G3 = -pS3./(pS3+2*pi*j*f);
h1=line_fewer_markers(f,20*log10(abs(G1)),8,'rs','spacing','logx','linewidth',2,'markerfacecolor','r')
h2=line_fewer_markers(f,20*log10(abs(G2)),8,'bp','spacing','logx','linewidth',2)
h3=line_fewer_markers(f,20*log10(abs(G3)),8,'mv','spacing','logx','linewidth',2,'markersize',6)
set(gca,'XScale','log')

legend([h1,h2,h3],'G1','G2','G3','location','southwest');
