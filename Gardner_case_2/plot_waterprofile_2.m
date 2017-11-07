%%
% clear; close all; clc
figure
x=linspace(0,-1,26);
f=fopen('g2_24_h.0000');h0=fread(f,26,'float64');fclose(f);
f=fopen('g2_24_h.0001');h1=fread(f,26,'float64');fclose(f);
f=fopen('g2_24_h.0002');h2=fread(f,26,'float64');fclose(f);
f=fopen('g2_24_h.0003');h3=fread(f,26,'float64');fclose(f);
f=fopen('g2_24_h.0004');h4=fread(f,26,'float64');fclose(f);
f=fopen('g2_24_h.0005');h5=fread(f,26,'float64');fclose(f);
f=fopen('g2_24_h.0006');h6=fread(f,26,'float64');fclose(f);
f=fopen('g2_24_h.0007');h7=fread(f,26,'float64');fclose(f);
plot(h0,x,'r','linewidth',1.5)
hold on
plot(h1,x)
plot(h2,x,'r','linewidth',1.5)
plot(h3,x,'r','linewidth',1.5)
plot(h4,x,'r','linewidth',1.5)
plot(h5,x,':r','linewidth',1.5)
plot(h6,x,'--r','linewidth',1.5)
plot(h7,x,'r','linewidth',2)

set(gca, 'fontsize',15,'ytick',[-1 -0.5 0],'yticklabel',[1 0.5 0])
grid on
xlabel('Pressure Head')
ylabel('Depth')
title('N=20','fontsize',12)


x=linspace(0,-1,76);
f=fopen('g2_74_h.0000');h0=fread(f,76,'float64');fclose(f);
f=fopen('g2_74_h.0001');h1=fread(f,76,'float64');fclose(f);
f=fopen('g2_74_h.0002');h2=fread(f,76,'float64');fclose(f);
f=fopen('g2_74_h.0003');h3=fread(f,76,'float64');fclose(f);
f=fopen('g2_74_h.0004');h4=fread(f,76,'float64');fclose(f);
f=fopen('g2_74_h.0005');h5=fread(f,76,'float64');fclose(f);
f=fopen('g2_74_h.0006');h6=fread(f,76,'float64');fclose(f);
f=fopen('g2_74_h.0007');h7=fread(f,76,'float64');fclose(f);
plot(h0,x,'k','linewidth',1.5)
hold on
plot(h1,x)
plot(h2,x,'k','linewidth',1.5)
plot(h3,x,'k','linewidth',1.5)
plot(h4,x,'k','linewidth',1.5)
plot(h5,x,':k','linewidth',1.5)
plot(h6,x,'--k','linewidth',1.5)
plot(h7,x,'k','linewidth',2)
set(gca, 'fontsize',15,'ytick',[-1 -0.5 0],'yticklabel',[1 0.5 0])
grid on
xlabel('Pressure Head')
ylabel('Depth')
title('N=74','fontsize',12)
