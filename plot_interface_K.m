%%
% close all
figure
% hold on
fdir=dir('g1_50_lin_time.bin');
n=fdir.bytes/8;
f=fopen(fdir.name);t50=fread(f,n,'float64');fclose(f);
f1=fopen('g1_50_lin_kplus.bin');kplus50=fread(f1,n,'float64');fclose(f1);
p1=loglog(t50,kminus50,':k','linewidth',1.5);
f2=fopen('g1_50_lin_kminus.bin');kminus50=fread(f2,n,'float64');fclose(f1);
hold on
loglog(t50,kplus50,':k','linewidth',1.5)

fdir=dir('g1_100_lin_time.bin');
n=fdir.bytes/8;
f=fopen(fdir.name);t100=fread(f,n,'float64');fclose(f);
f1=fopen('g1_100_lin_kplus.bin');kplus100=fread(f1,n,'float64');fclose(f1);
p2=loglog(t100,kminus100,'--k','linewidth',1.5);
f2=fopen('g1_100_lin_kminus.bin');kminus100=fread(f2,n,'float64');fclose(f1);
loglog(t100,kplus100,'--k','linewidth',1.5)



fdir=dir('g1_200_lin_time.bin');
n=fdir.bytes/8;
f=fopen(fdir.name);t200=fread(f,n,'float64');fclose(f);
f1=fopen('g1_200_lin_kplus.bin');kplus200=fread(f1,n,'float64');fclose(f1);
p3=loglog(t200,kminus200,'-k','linewidth',1);
f2=fopen('g1_200_lin_kminus.bin');kminus200=fread(f2,n,'float64');fclose(f1);
loglog(t200,kplus200,'-k','linewidth',1)

% title('Gardner Model')
set(gca,'fontsize',15)
xlabel('Time','interpreter','latex')
ylabel('Conductivity at the interface','interpreter','latex')
le=legend([p1 p2 p3],'$$N=50$','$N=100$',...
    '$N=200$',...
    'location','best');
set(le,'interpreter','latex');

% xlim([0, 10])
xlim([1e-1, 1e2])
ylim([1e-8, 1e-3])

text(1,10^(-5.3),'$K^+$','interpreter','latex','fontsize',15)

text(1,10^(-3.5),'$K^-$','interpreter','latex','fontsize',15)
