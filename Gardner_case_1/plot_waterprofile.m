%%
clear; close all; clc
figure
% subplot(2,2,1)
x=linspace(0,-1,52);
f=fopen('g1_50_lin_h.0000');h0=fread(f,52,'float64');fclose(f);
% f=fopen('g1_50_lin_h.0001');h1=fread(f,52,'float64');fclose(f);
f=fopen('g1_50_lin_h.0002');h2=fread(f,52,'float64');fclose(f);
f=fopen('g1_50_lin_h.0003');h3=fread(f,52,'float64');fclose(f);
f=fopen('g1_50_lin_h.0004');h4=fread(f,52,'float64');fclose(f);
f=fopen('g1_50_lin_h.0005');h5=fread(f,52,'float64');fclose(f);
f=fopen('g1_50_lin_h.0006');h6=fread(f,52,'float64');fclose(f);
f=fopen('g1_50_lin_h.0007');h7=fread(f,52,'float64');fclose(f);
plot(h0,x,'k','linewidth',1.5)
hold on
% plot(h1,x)
plot(h2,x,'k','linewidth',1.5)
plot(h3,x,'k','linewidth',1.5)
plot(h4,x,'k','linewidth',1.5)
plot(h5,x,':k','linewidth',1.5)
plot(h6,x,'--k','linewidth',1.5)
plot(h7,x,'k','linewidth',2)
xlim([-1.65 0])
text(-1.3,-0.04,'$t=0$','interpreter','latex','fontsize',14)
text(-1.2,-0.11,'$t=0.001$','interpreter','latex','fontsize',14)
text(-1.05,-0.26,'$t=0.01$','interpreter','latex','fontsize',14)
text(-0.93,-0.32,'$t=0.05$','interpreter','latex','fontsize',14)
text(-0.65,-0.3,'$t=2.5$','interpreter','latex','fontsize',14)
text(-0.68,-0.35,'$t=4.5$','interpreter','latex','fontsize',14);
text(-0.22,-0.34,'$t=100$','interpreter','latex','fontsize',14);
text(-1.44,-0.7,'$t=100$','interpreter','latex','fontsize',14);
txx = text(-0.53,-0.35,'$\longrightarrow$','interpreter','latex','fontsize',20);
set(txx,'rotation',25)
set(gca, 'fontsize',15,'ytick',[-1 -0.5 0],'yticklabel',[1 0.5 0])
grid on
xlabel('Pressure Head','interpreter','latex')
ylabel('Depth','interpreter','latex')
% title('$N=50$ with multiple solutions (permanently) to the interface problem','interpreter','latex','fontsize',14)
%%
% subplot(2,2,2)
x=linspace(0,-1,102);
f=fopen('g1_100_lin_h.0000');h0=fread(f,102,'float64');fclose(f);
% f=fopen('g1_100_lin_h.0001');h1=fread(f,102,'float64');fclose(f);
f=fopen('g1_100_lin_h.0002');h2=fread(f,102,'float64');fclose(f);
f=fopen('g1_100_lin_h.0003');h3=fread(f,102,'float64');fclose(f);
f=fopen('g1_100_lin_h.0004');h4=fread(f,102,'float64');fclose(f);
f=fopen('g1_100_lin_h.0005');h5=fread(f,102,'float64');fclose(f);
f=fopen('g1_100_lin_h.0006');h6=fread(f,102,'float64');fclose(f);
f=fopen('g1_100_lin_h.0007');h7=fread(f,102,'float64');fclose(f);
plot(h0,x,'k','linewidth',1.5)
hold on
% plot(h1,x)
plot(h2,x,'k','linewidth',1.5)
plot(h3,x,'k','linewidth',1.5)
plot(h4,x,'k','linewidth',1.5)
plot(h5,x,':k','linewidth',1.5)
plot(h6,x,'--k','linewidth',1.5)
plot(h7,x,'k','linewidth',2)
xlim([-1.65 0])
text(-1.4,-0.04,'$t=0$','interpreter','latex','fontsize',14)
text(-1.2,-0.12,'$t=0.001$','interpreter','latex','fontsize',14)
text(-1.1,-0.28,'$t=0.01$','interpreter','latex','fontsize',14)
text(-0.98,-0.35,'$t=0.05$','interpreter','latex','fontsize',14)
text(-0.4,-0.4,'$t=2.5$','interpreter','latex','fontsize',14)
text(-0.57,-0.45,'$t=4.5$','interpreter','latex','fontsize',14);
% text(-0.766,-0.39,'$t=100$','interpreter','latex','fontsize',14);
text(-0.8,-0.7,'$t=100$','interpreter','latex','fontsize',14);
set(gca, 'fontsize',15,'ytick',[-1 -0.5 0],'yticklabel',[1 0.5 0])
grid on
xlabel('Pressure Head','interpreter','latex')
ylabel('Depth','interpreter','latex')
% title('$N=100$ with multiple solutions (temporarily) to the interface problem','interpreter','latex','fontsize',14)

%%
% subplot(2,2,3)
x=linspace(0,-1,202);
f=fopen('g1_200_lin_h.0000');h0=fread(f,202,'float64');fclose(f);
% f=fopen('g1_200_lin_h.0001');h1=fread(f,202,'float64');fclose(f);
f=fopen('g1_200_lin_h.0002');h2=fread(f,202,'float64');fclose(f);
f=fopen('g1_200_lin_h.0003');h3=fread(f,202,'float64');fclose(f);
f=fopen('g1_200_lin_h.0004');h4=fread(f,202,'float64');fclose(f);
f=fopen('g1_200_lin_h.0005');h5=fread(f,202,'float64');fclose(f);
f=fopen('g1_200_lin_h.0006');h6=fread(f,202,'float64');fclose(f);
f=fopen('g1_200_lin_h.0007');h7=fread(f,202,'float64');fclose(f);
plot(h0,x,'k','linewidth',1.5)
hold on
% plot(h1,x)
plot(h2,x,'k','linewidth',1.5)
plot(h3,x,'k','linewidth',1.5)
plot(h4,x,'k','linewidth',1.5)
plot(h5,x,':k','linewidth',1.5)
plot(h6,x,'--k','linewidth',1.5)
plot(h7,x,'k','linewidth',2)
xlim([-1.65 0])
text(-1.3,-0.03,'$t=0$','interpreter','latex','fontsize',14)
text(-1.2,-0.12,'$t=0.001$','interpreter','latex','fontsize',14)
text(-1.1,-0.28,'$t=0.01$','interpreter','latex','fontsize',14)
text(-0.95,-0.33,'$t=0.05$','interpreter','latex','fontsize',14)
text(-0.88,-0.48,'$t=2.5$','interpreter','latex','fontsize',14)
text(-0.96,-0.58,'$t=4.5$','interpreter','latex','fontsize',14);
text(-0.8,-0.7,'$t=100$','interpreter','latex','fontsize',14)
set(gca, 'fontsize',15,'ytick',[-1 -0.5 0],'yticklabel',[1 0.5 0])
grid on
xlabel('Pressure Head','interpreter','latex')
ylabel('Depth','interpreter','latex')
% title('$N=200$ with a unique solution to the interface problem','interpreter','latex','fontsize',14)

%%
% subplot(2,2,4)
n=6016;
f=fopen('g1_50_lin_time.bin');t05=fread(f,n,'float64');fclose(f);
f=fopen('g1_50_lin_hbelow_neq50.bin');hb05=fread(f,n,'float64');fclose(f);
semilogx(t05,hb05,':k','linewidth',1.5)
n=6058;
f=fopen('g1_100_lin_time.bin');t1=fread(f,n,'float64');fclose(f);
f=fopen('g1_100_lin_hbelow_neq50.bin');hb1=fread(f,n,'float64');fclose(f);
hold on, semilogx(t1,hb1,'--k','linewidth',1.5)
n=6021;
f=fopen('g1_200_lin_time.bin');t2=fread(f,n,'float64');fclose(f);
f=fopen('g1_200_lin_hbelow_neq50.bin');hb2=fread(f,n,'float64');fclose(f);
hold on, semilogx(t2,hb2,'k','linewidth',1.5)

set(gca,'fontsize',15)
xlabel('Time','interpreter','latex')
ylabel('Pressure Head at $z \approx 0.5098$','interpreter','latex')
xlim([1e-1 1e2])
text(1.7,-1.3,'$N=50$ with multiple solutions $\rightarrow$','interpreter','latex','fontsize',14)
text(2.8,-1.0,'$\leftarrow$ $N=100$ with multiple solutions','interpreter','latex','fontsize',14)
text(1.3,-0.8,'$N=200$ with a unique solution','interpreter','latex','fontsize',14)
txx = text(28,-0.8,'$\rightarrow$','interpreter','latex','fontsize',14);
set(txx,'rotation',-30)

n=2728;
f=fopen('g1_50_lin_multitime.bin');mt=fread(f,n,'float64');fclose(f);
y=zeros(1,1);
for i=1
    y(i)=(hb05(t05==mt(i)));
end
hold on, semilogx(mt(1),y,'ok')
n=1424;
f=fopen('g1_100_lin_multitime.bin');mt=fread(f,n,'float64');fclose(f);
y=zeros(1,1);
for i=1
    y(i)=(hb1(t1==mt(i)));
end
hold on, semilogx(mt(1),y,'ok')
y=zeros(1,1);
for i=n
    y(i)=(hb1(t1==mt(i)));
end
hold on, semilogx(mt(n),y,'ok')

tx = text(0.73,-1.04,'\leftarrow','fontsize',15);
set(tx,'rotation',-90)
text(0.2,-1.09,'multiple solutions appear here','interpreter','latex','fontsize',14);
tx2 = text(1.03,-1.055,'$\rightarrow$','interpreter','latex','fontsize',14);
set(tx2,'rotation',90)
tx3 = text(3.4,-0.895,'$\rightarrow$','interpreter','latex','fontsize',14);
set(tx3,'rotation',-50)
text(0.2,-0.88,'multiple solutions disappear here','interpreter','latex','fontsize',14);

ylim([-1.4,-0.7])
ld=legend('$N=50$','$N=100$','$N=200$','location','southwest');
set(ld,'interpreter','latex');
shg
