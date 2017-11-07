%%
clear; close all; clc
figure
x=linspace(0,-1,402);
f=fopen('f_s_4_r2_h.0000');h0=fread(f,402,'float64');fclose(f);
% f=fopen('f_s_4_r2_h.0001');h1=fread(f,402,'float64');fclose(f);
f=fopen('f_s_4_r2_h.0002');h2=fread(f,402,'float64');fclose(f);
f=fopen('f_s_4_r2_h.0003');h3=fread(f,402,'float64');fclose(f);
f=fopen('f_s_4_r2_h.0004');h4=fread(f,402,'float64');fclose(f);
% f=fopen('f_s_4_r2_h.0005');h5=fread(f,402,'float64');fclose(f);
% f=fopen('f_s_4_r2_h.0006');h6=fread(f,402,'float64');fclose(f);
f=fopen('f_s_4_r2_h.0007');h7=fread(f,402,'float64');fclose(f);
p1 = plot(h0,x,':k','linewidth',2)
hold on
% plot(h1,x)
plot(h2,x,':k','linewidth',2)
plot(h3,x,':k','linewidth',2)
plot(h4,x,':k','linewidth',2)
% plot(h5,x,'--k','linewidth',1.5)
% plot(h6,x,'k','linewidth',1.5)
plot(h7,x,':k','linewidth',2)
set(gca, 'fontsize',15,'ytick',[-1 -0.5 0],'yticklabel',[1 0.5 0])
grid on
xlabel('Pressure Head')
ylabel('Depth')
% title('N=400','fontsize',12)
ylim([-1.1 0])
xlim([-75 -35])
%%
x=linspace(0,-1,802);
f=fopen('f_s_8_r2_h.0000');h0=fread(f,802,'float64');fclose(f);
f=fopen('f_s_8_r2_h.0001');h1=fread(f,802,'float64');fclose(f);
f=fopen('f_s_8_r2_h.0002');h2=fread(f,802,'float64');fclose(f);
f=fopen('f_s_8_r2_h.0003');h3=fread(f,802,'float64');fclose(f);
f=fopen('f_s_8_r2_h.0004');h4=fread(f,802,'float64');fclose(f);
% f=fopen('f_s_8_r2_h.0005');h5=fread(f,802,'float64');fclose(f);
% f=fopen('f_s_8_r2_h.0006');h6=fread(f,802,'float64');fclose(f);
f=fopen('f_s_8_r2_h.0007');h7=fread(f,802,'float64');fclose(f);
p2 = plot(h0,x,'--k','linewidth',1.5)
hold on
% plot(h1,x)
plot(h2,x,'--k','linewidth',1.5)
plot(h3,x,'--k','linewidth',1.5)
plot(h4,x,'--k','linewidth',1.5)
% plot(h5,x,'--k','linewidth',1.5)
% plot(h6,x,'k','linewidth',1.5)
plot(h7,x,'--k','linewidth',2)
text(-63,-0.04,'t=0','fontsize',15)
text(-63,-0.17,'t=10','fontsize',15)
text(-63,-0.51,'t=50','fontsize',15)
text(-63,-0.8,'t=100','fontsize',15)
% text(-46,-0.9,'t=150','fontsize',15)
text(-63,-1.03,'t=1000','fontsize',15);
set(gca, 'fontsize',15,'ytick',[-1 -0.5 0],'yticklabel',[1 0.5 0])
grid on
xlabel('Pressure Head')
ylabel('Depth')
% title('N=800','fontsize',12)

% %%
% x=linspace(0,-1,1602);
% f=fopen('f_s_16_r2_h.0000');h0=fread(f,1602,'float64');fclose(f);
% f=fopen('f_s_16_r2_h.0001');h1=fread(f,1602,'float64');fclose(f);
% f=fopen('f_s_16_r2_h.0002');h2=fread(f,1602,'float64');fclose(f);
% f=fopen('f_s_16_r2_h.0003');h3=fread(f,1602,'float64');fclose(f);
% f=fopen('f_s_16_r2_h.0004');h4=fread(f,1602,'float64');fclose(f);
% f=fopen('f_s_16_r2_h.0005');h5=fread(f,1602,'float64');fclose(f);
% f=fopen('f_s_16_r2_h.0006');h6=fread(f,1602,'float64');fclose(f);
% f=fopen('f_s_16_r2_h.0007');h7=fread(f,1602,'float64');fclose(f);
% plot(h0,x,'k','linewidth',1.5)
% hold on
% % plot(h1,x)
% plot(h2,x,'k','linewidth',1.5)
% plot(h3,x,'k','linewidth',1.5)
% plot(h4,x,'k','linewidth',1.5)
% % plot(h5,x,'k','linewidth',1.5)
% % plot(h6,x,'--k','linewidth',1.5)
% plot(h7,x,'k','linewidth',2)
% set(gca, 'fontsize',15,'ytick',[-1 -0.5 0],'yticklabel',[1 0.5 0])
% grid on
% xlabel('Pressure Head')
% ylabel('Depth')
% % title('N=1600','fontsize',12)

%%
x=linspace(0,-1,3202);
f=fopen('f_s_32_r2_h.0000');h0=fread(f,3202,'float64');fclose(f);
f=fopen('f_s_32_r2_h.0001');h1=fread(f,3202,'float64');fclose(f);
f=fopen('f_s_32_r2_h.0002');h2=fread(f,3202,'float64');fclose(f);
f=fopen('f_s_32_r2_h.0003');h3=fread(f,3202,'float64');fclose(f);
f=fopen('f_s_32_r2_h.0004');h4=fread(f,3202,'float64');fclose(f);
f=fopen('f_s_32_r2_h.0005');h5=fread(f,3202,'float64');fclose(f);
f=fopen('f_s_32_r2_h.0006');h6=fread(f,3202,'float64');fclose(f);
f=fopen('f_s_32_r2_h.0007');h7=fread(f,3202,'float64');fclose(f);
p3 = plot(h0,x,'k','linewidth',1.5)
hold on
% plot(h1,x)
plot(h2,x,'k','linewidth',1.5)
plot(h3,x,'k','linewidth',1.5)
plot(h4,x,'k','linewidth',1.5)
% plot(h5,x,'k','linewidth',1.5)
% plot(h6,x,'--k','linewidth',1.5)
plot(h7,x,'k','linewidth',2)
set(gca, 'fontsize',15,'ytick',[-1 -0.5 0],'yticklabel',[1 0.5 0])
grid on
xlabel('Pressure Head','interpreter','latex')
ylabel('Depth','interpreter','latex')
% title('N=1600','fontsize',12)

% legend([p1,p2,p3],'N=400, multiple solutions',...
%     'N=800, multiple solutions','fine mesh N=3200, unique solution',...
%     'location','southoutside','Orientation','vertical')

le = legend([p1,p2,p3],'$N=400$',...
    '$N=800$','$N=3200$',...
    'location','best');
set(le,'interpreter','latex')