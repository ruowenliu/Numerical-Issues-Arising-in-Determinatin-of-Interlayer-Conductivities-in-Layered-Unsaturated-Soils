%% Gardner
close all;clc;
param1.a=13;param2.a=1;param1.b=1;param2.b=6e-4;param1.tr=0.06;param2.tr=0.06;param1.ts=0.4;param2.ts=0.4;
x=logspace(-3,2,101);


sat1=ShG(-x,param1);sat2=ShG(-x,param2);
semilogx(x,sat1,'--k',x,sat2,'+k')
set(gca, 'FontSize', 15)
xlabel('Soil Suction (L)','interpreter','latex')
ylabel('Effective Saturation','interpreter','latex')
ld1=legend('upper layer','lower layer','location','northeast');
set(ld1,'interpreter','latex')
set(gca,'xtick',[1e-3 1e-2 1e-1 1 1e1 1e2])
ylim([0,1.1])
grid on

figure
k1=KSG(sat1,param1);k2=KSG(sat2,param2);
loglog(x,k1,'--k',x,k2,'+k')
set(gca, 'FontSize', 15)
xlabel('Soil Suction (L)','interpreter','latex')
ylabel('Hydraulic Conductivity (L)','interpreter','latex')
ld2=legend('upper layer','lower layer','location','southwest');
set(ld2,'interpreter','latex')
set(gca,'xtick',[1e-3 1e-2 1e-1 1 1e1 1e2])
ylim([1e-40 10])
grid on

