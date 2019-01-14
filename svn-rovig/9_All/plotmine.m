figure
hold on

load ('squaren1C1F2Lambda50.mat')
resSquare{1}=residual;
load ('squaren1C1F3Lambda50.mat')
resSquare{2}=residual;
load ('squaren1C1F4Lambda50.mat')
resSquare{3}=residual;
load ('squaren1C1F5Lambda50.mat')
resSquare{4}=residual;
load ('squaren1C1F6Lambda50.mat')
resSquare{5}=residual;
load ('squaren1C1F7Lambda50.mat')
resSquare{6}=residual;


for lev=1:6
minimum(lev)=min(resSquare{lev});
end
minimum=max(minimum);
minimum=10^(-11);
for lev=1:6
  resSquare{lev}(resSquare{lev}<=minimum)=[];
end
maximum=25;
% plot(log10(resSquare{1}(1:min(length(resSquare{1}),maximum))),'b-','LineWidth',12);
% plot(log10(resSquare{2}(1:min(length(resSquare{2}),maximum))),'g-','LineWidth',12);
% plot(log10(resSquare{3}(1:min(length(resSquare{3}),maximum))),'r-','LineWidth',12);
% plot(log10(resSquare{4}(1:min(length(resSquare{4}),maximum))),'k-','LineWidth',12);
% plot(log10(resSquare{5}(1:min(length(resSquare{5}),maximum))),'c-','LineWidth',12);
% plot(log10(resSquare{6}(1:min(length(resSquare{6}),maximum))),'m-','LineWidth',12);



plot(resSquare{1}(2:min(length(resSquare{1}),maximum))./resSquare{1}(1:min(length(resSquare{1}),maximum)-1),'b-','LineWidth',12);
plot(resSquare{2}(2:min(length(resSquare{2}),maximum))./resSquare{2}(1:min(length(resSquare{2}),maximum)-1),'g-','LineWidth',12);
plot(resSquare{3}(2:min(length(resSquare{3}),maximum))./resSquare{3}(1:min(length(resSquare{3}),maximum)-1),'r-','LineWidth',12);
plot(resSquare{4}(2:min(length(resSquare{4}),maximum))./resSquare{4}(1:min(length(resSquare{4}),maximum)-1),'k-','LineWidth',12);
plot(resSquare{5}(2:min(length(resSquare{5}),maximum))./resSquare{5}(1:min(length(resSquare{5}),maximum)-1),'c-','LineWidth',12);
plot(resSquare{6}(2:min(length(resSquare{6}),maximum))./resSquare{6}(1:min(length(resSquare{6}),maximum)-1),'m-','LineWidth',12);


l=legend('L=2, Dofs=50 ','L=3, Dofs=162 ','L=4, Dofs=578 ','L=5, Dofs=2178 ', 'L=6, Dofs=8450', 'L=7, Dofs=33282' );
l.FontSize=46;
xx=xlabel('Iterations');
yy=ylabel('log10(Residual)');
set(gca,'fontsize',48)
title('\lambda=1e50, \mu=1')

ylim([-9.5 -4])






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
