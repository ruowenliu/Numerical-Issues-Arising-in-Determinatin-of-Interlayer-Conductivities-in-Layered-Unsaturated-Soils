%%
figure
n=6843;
f=fopen('g2_24_time.bin');t01=fread(f,n,'float64');fclose(f);
f=fopen('g2_24_hbelow_neq24.bin');hb01=fread(f,n,'float64');fclose(f);
semilogx(t01,hb01,'--k','linewidth',1.5)
hold on

n=89282;
f=fopen('g2_24_200_time.bin');t02=fread(f,n,'float64');fclose(f);
f=fopen('g2_24_200_hbelow_neq24.bin');hb02=fread(f,n,'float64');fclose(f);
semilogx(t02,hb02,':k','linewidth',1.8)
set(gca,'fontsize',15)


n=6551;
f=fopen('g2_24_nt_time.bin');t06=fread(f,n,'float64');fclose(f);
f=fopen('g2_24_nt_hbelow_neq24.bin');hb06=fread(f,n,'float64');fclose(f);
hold on, loglog(t06,hb06,'k','linewidth',1)

ld = legend('Fixed-point iteration, \textit{maxiter} $ =199$',...
    'Fixed-point iteration, \textit{maxiter} $ =200$',...
    'Newton iteration','location','southwest');
set(ld,'interpreter','latex')


n=179;
f=fopen('g2_24_fxptdivtime.bin');dvt=fread(f,n,'float64');fclose(f);
y=(hb01(t01==dvt(1)));
hold on, semilogx(dvt(1),y,'.k','MarkerSize',30)
y=(hb01(t01==dvt(end)));
hold on, semilogx(dvt(end),y,'.k','MarkerSize',30)

n=82867;
f=fopen('g2_24_200_fxptdivtime.bin');dvt=fread(f,n,'float64');fclose(f);
y=(hb02(t02==dvt(1)));
hold on, semilogx(dvt(1),y,'.k','MarkerSize',30)
y=(hb02(t02==dvt(end)));
hold on, semilogx(dvt(end),y,'.k','MarkerSize',30)

% n=6456;
% f=fopen('g2_74_time.bin');t06=fread(f,n,'float64');fclose(f);
% f=fopen('g2_74_hbelow_neq24.bin');hb06=fread(f,n,'float64');fclose(f);
% hold on, loglog(t06,hb06,':k','linewidth',1.5)

% n=6532;
% f=fopen('g2_224_time.bin');t06=fread(f,n,'float64');fclose(f);
% f=fopen('g2_224_hbelow_neq24.bin');hb06=fread(f,n,'float64');fclose(f);
% hold on, loglog(t06,hb06,':k','linewidth',1.5)

set(gca,'fontsize',15)
xlabel('Time','interpreter','latex')
ylabel('Pressure Head at $z \approx 0.52$','interpreter','latex')
xlim([1e-8 1e2])
ylim([-.12 -0.01])
text(3.2e-5,-0.0457,'$\leftarrow$','interpreter','latex','fontsize',25)
text(1.5e-4,-0.0457,'fixed-point: 2-cycle appears','interpreter','latex','fontsize',15)
text(0.002,-0.056,'2-cycle disappears','interpreter','latex','fontsize',15)
tx1=text(0.006,-0.065,'$\leftarrow$','interpreter','latex','fontsize',25);
set(tx1,'rotation',45)
tx2=text(0.0032,-0.065,'$\leftarrow$','interpreter','latex','fontsize',25);
set(tx2,'rotation',45)

% xlim([2e-5 1e-2])
% ylim([-.07 -0.04])