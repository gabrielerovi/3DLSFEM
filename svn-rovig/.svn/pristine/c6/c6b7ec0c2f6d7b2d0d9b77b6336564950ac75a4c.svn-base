figure
hold on

load ('circleresn45m25L2TBminMidpoint.mat')
resn45m25L2TBminMidpoint=residual;
load ('circleresn35m25L2TBminMidpoint.mat')
resn35m25L2TBminMidpoint=residual;

load ('circleresn45m25L2TBmin.mat')
resn45m25L2TBmin=residual;
load ('circleresn35m25L2TBmin.mat')
resn35m25L2TBmin=residual;

load ('circleresn45m25L2TB.mat')
resn45m25L2TB=residual;
load ('circleresn35m25L2TB.mat')
resn35m25L2TB=residual;


load ('circleresZn45m25L2TB.mat')
resZn45m25L2TB=resz(1:50);
load('circleresZn35m25L2TB.mat');
resZn35m25L2TB=resz(1:50);

plot(log10(resn35m25L2TBmin),'bd','MarkerSize',20,'MarkerFaceColor','blue');
plot(log10(resn45m25L2TBmin),'b-','LineWidth',12);
plot(log10(resn35m25L2TBminMidpoint),'kd','MarkerSize',20,'MarkerFaceColor','black');
plot(log10(resn45m25L2TBminMidpoint),'k-','LineWidth',12);
plot(log10(resn35m25L2TB),'rd','MarkerSize',20,'MarkerFaceColor','red');
plot(log10(resn45m25L2TB),'r-','LineWidth',12);

plot(log10(resZn45m25L2TB(1:50)),'g-','LineWidth',6);
l=legend('I) min','II) min','I) midpoint','II) midpoint','I) linear','II) linear');
l.FontSize=46;
xx=xlabel('Iterations');
yy=ylabel('log10(Residual)');
set(gca,'fontsize',48)



figure 

hold on

load('circleresn25m25L23TB.mat');
resn25m25L23TB=residual;
load('circleresZn25m25L23TB.mat');
resZn25m25L23TB=resz(1:50);

load('circleresn25m25L23TBmin.mat');
resn25m25L23TBmin=residual;
load('circleresZn25m25L23TBmin.mat');
resZn25m25L23TBmin=resz(1:50);

load('circleresn25m25L13TB.mat');
resn25m25L13TB=residual;
load('circleresZn25m25L13TB.mat');
resZn25m25L13TB=resz(1:50);


load('circleresn25m25L13TBmin.mat');
resn25m25L13TBmin=residual;
load('circleresZn25m25L13TBmin.mat');
resZn25m25L13TBmin=resz(1:50);


plot(log10(resn25m25L13TBmin),'bd','MarkerSize',20,'MarkerFaceColor','blue');
plot(log10(resn25m25L23TBmin),'b-','LineWidth',12);
plot(log10( resn25m25L13TB),'rd','MarkerSize',20,'MarkerFaceColor','red');
plot(log10( resn25m25L23TB),'r-','LineWidth',12);


plot(log10(resZn25m25L13TB),'g-','LineWidth',6);
plot(log10(resZn25m25L23TB),'g-','LineWidth',6);

l=legend('3L min','2L min','3L linear','2L linear');
l.FontSize=46;
xx=xlabel('Iterations');
yy=ylabel('log10(Residual)');
set(gca,'fontsize',48)

