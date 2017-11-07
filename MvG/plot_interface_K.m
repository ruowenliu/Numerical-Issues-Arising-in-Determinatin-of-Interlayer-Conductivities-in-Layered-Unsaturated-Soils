%%
% close all
figure
fdir=dir('m_104_r2_time.bin');
n=fdir.bytes/8;
f=fopen(fdir.name);t104=fread(f,n,'float64');fclose(f);
f1=fopen('m_104_r2_kplus.bin');kplus104=fread(f1,n,'float64');fclose(f1);
p1=semilogy(t104,kplus104,'--k','linewidth',1.5);
f2=fopen('m_104_r2_kminus.bin');kminus104=fread(f2,n,'float64');fclose(f1);
hold on
semilogy(t104,kminus104,'--k','linewidth',1.5)


fdir=dir('m_110_r0_time.bin');
n=fdir.bytes/8;
f=fopen(fdir.name);t110=fread(f,n,'float64');fclose(f);
f1=fopen('m_110_r0_kplus.bin');kplus110=fread(f1,n,'float64');fclose(f1);
p2=semilogy(t110,kplus110,'-k','linewidth',1.5);
f2=fopen('m_110_r0_kminus.bin');kminus110=fread(f2,n,'float64');fclose(f1);
hold on
semilogy(t110,kminus110,'-k','linewidth',1.5)

%%

set(gca,'fontsize',15)
xlabel('Time','interpreter','latex')
ylabel('Conductivity at the interface','interpreter','latex')
le=legend([p1 p2],'$N=104$','$N=110$',...
    'location','best');
set(le,'interpreter','latex');

% xlim([1e-1, 1e2])
xlim([0.014, 0.016])
text(0.0142,10^(-10.2),'$K^+$','interpreter','latex','fontsize',15)
text(0.0142,10^(-7.5),'$K^-$','interpreter','latex','fontsize',15)
