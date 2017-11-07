%%
% close all
rc=3; 

figure
fdir=dir(['f_s_4_r' num2str(rc) '_time.bin']);
n=fdir.bytes/8;
f=fopen(fdir.name);t400=fread(f,n,'float64');fclose(f);
f1=fopen(['f_s_4_r' num2str(rc) '_kplus.bin']);kplus400=fread(f1,n,'float64');fclose(f1);
p1=semilogy(t400,kplus400,':k','linewidth',1.5);
f2=fopen(['f_s_4_r' num2str(rc) '_kminus.bin']);kminus400=fread(f2,n,'float64');fclose(f1);
hold on
% hold off
% p1=semilogy(t400,kminus400,':k','linewidth',1.5);


fdir=dir(['f_s_8_r' num2str(rc) '_time.bin']);
n=fdir.bytes/8;
f=fopen(fdir.name);t800=fread(f,n,'float64');fclose(f);
f1=fopen(['f_s_8_r' num2str(rc) '_kplus.bin']);kplus800=fread(f1,n,'float64');fclose(f1);
hold on
p2=semilogy(t800,kplus800,'--k','linewidth',1.5);
f2=fopen(['f_s_8_r' num2str(rc) '_kminus.bin']);kminus800=fread(f2,n,'float64');fclose(f1);
hold on
% semilogy(t800,kminus800,'--k','linewidth',1.5)

rc=2;
fdir=dir(['f_s_32_r' num2str(rc) '_time.bin']);
n=fdir.bytes/8;
f=fopen(fdir.name);t1600=fread(f,n,'float64');fclose(f);
f1=fopen(['f_s_32_r' num2str(rc) '_kplus.bin']);kplus1600=fread(f1,n,'float64');fclose(f1);
p3=semilogy(t1600,kplus1600,'-k','linewidth',1.5);
f2=fopen(['f_s_32_r' num2str(rc) '_kminus.bin']);kminus1600=fread(f2,n,'float64');fclose(f1);
hold on
% semilogy(t1600,kminus1600,'-k','linewidth',1.5)

%%

set(gca,'fontsize',15)
xlabel('Time','interpreter','latex')
ylabel('Conductivity at the interface','interpreter','latex')
le=legend([p1 p2 p3],'$N=400$','$N=800$','$N=3200$',...
    'location','best');
set(le,'interpreter','latex');
% switch rc
%     case 0
%         title(['FXLR Model, $r_0 =$ last root'], 'interpreter','latex')
%     case 1
%         title(['FXLR Model, $r_0 =$ leftmost'], 'interpreter','latex')
%     case 2
%         title(['FXLR Model, $r_0 =$ middle'], 'interpreter','latex')
%     case 3
%         title(['FXLR Model, $r_0 =$ rightmost'], 'interpreter','latex')
% end

xlim([30, 300])
ylim([10^(-20), 10])
% text(0.0142,10^(-10.2),'$K^+$','interpreter','latex','fontsize',15)
% text(0.0142,10^(-7.5),'$K^-$','interpreter','latex','fontsize',15)

