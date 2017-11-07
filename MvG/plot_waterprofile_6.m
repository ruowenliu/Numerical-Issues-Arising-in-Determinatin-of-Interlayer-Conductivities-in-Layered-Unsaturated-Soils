close all; clear; clc
figure
subplot(2,1,2)
% n=110 unique root
n=110;
x=linspace(0,-1,n+2);
% f=fopen(['m_' num2str(n) '_loc_w.0000']);w0=fread(f,n+2,'float64');fclose(f);
% f=fopen(['m_' num2str(n) '_loc_w.0001']);w1=fread(f,n+2,'float64');fclose(f);
% f=fopen(['m_' num2str(n) '_loc_w.0002']);w2=fread(f,n+2,'float64');fclose(f);
f=fopen(['m_' num2str(n) '_r0_w.0003']);w3=fread(f,n+2,'float64');fclose(f);
f=fopen(['m_' num2str(n) '_r0_w.0004']);w4=fread(f,n+2,'float64');fclose(f);
f=fopen(['m_' num2str(n) '_r0_w.0005']);w5=fread(f,n+2,'float64');fclose(f);
f=fopen(['m_' num2str(n) '_r0_w.0006']);w6=fread(f,n+2,'float64');fclose(f);

% p0 = plot(w0,x,'-k','linewidth',2);
hold on
% p1 = plot(w1,x,'-k','linewidth',2);
% p2 = plot(w2,x,'-k','linewidth',2);
p3 = plot(w3,x,'-k','linewidth',2);
% p4 = plot(w4,x,'.-k','linewidth',1);
p5 = plot(w5,x,'--k','linewidth',2);
p6 = plot(w6,x,':k','linewidth',2);
l = legend('$t=0.0149$','$t=0.0151$','$t=0.0152$','location','southeast');
set(l,'interpreter','latex')

set(gca, 'fontsize',15,'ytick',[-1 -0.5 0],'yticklabel',[1 0.5 0])
grid on
xlabel('Water Content','interpreter','latex')
ylabel('Depth','interpreter','latex')
% title('N=100','fontsize',12)
% legend([p0,p1,p2,p3,p4],'t=0','t=0.015','t=0.0151','t=0.0152','t=0.05');
title('$N=110$, unique solution','interpreter','latex')
xlim([0,0.25])
ylim([-0.515 -0.47])

%%
subplot(2,1,1)
n=104; % multiple roots
x=linspace(0,-1,n+2);
% f=fopen(['m_' num2str(n) '_loc_w.0000']);w0=fread(f,n+2,'float64');fclose(f);
% f=fopen(['m_' num2str(n) '_loc_w.0001']);w1=fread(f,n+2,'float64');fclose(f);
% f=fopen(['m_' num2str(n) '_loc_w.0002']);w2=fread(f,n+2,'float64');fclose(f);
f=fopen(['m_' num2str(n) '_r2_w.0003']);w3=fread(f,n+2,'float64');fclose(f);
f=fopen(['m_' num2str(n) '_r2_w.0004']);w4=fread(f,n+2,'float64');fclose(f);
f=fopen(['m_' num2str(n) '_r2_w.0005']);w5=fread(f,n+2,'float64');fclose(f);
f=fopen(['m_' num2str(n) '_r2_w.0006']);w6=fread(f,n+2,'float64');fclose(f);

% p00 = plot(w0,x,':k','linewidth',2);
hold on
% p11 = plot(w1,x,':k','linewidth',2);
% p22 = plot(w2,x,':k','linewidth',2);
p33 = plot(w3,x,'-k','linewidth',2);
% p44 = plot(w4,x,'.k','linewidth',2);
p55 = plot(w5,x,'--k','linewidth',2);

p66 = plot(w6,x,':k','linewidth',2);
% l=legend('$t=0.0149$','$t=0.0151$','$t=0.0152$','location','southeast');
% set(l,'interpreter','latex')

set(gca, 'fontsize',15,'ytick',[-1 -0.5 0],'yticklabel',[1 0.5 0])
grid on

ylabel('Depth','interpreter','latex')
title('$N=104$, multiple solutions','interpreter','latex')

xlim([0,0.25])
ylim([-0.515 -0.47])

% text(0.1,-0.04,'t=0','fontsize',15)
% text(0.1,-0.27,'t=0.007','fontsize',15)
% text(0.1,-0.44,'t=0.014','fontsize',15)
% text(0.1,-0.83,'t=0.05','fontsize',15)

% legend([p00,p0],'N=104, multiple solutions','N=110, unique solution','location','southoutside')

% title('N=150','fontsize',12)
% legend('t=0','t=0.014','t=0.015','t=0.0151','t=0.0152','t=0.05','t=0','t=0.014','t=0.015','t=0.0151','t=0.0152','t=0.05');
