%%
close all
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