close all
clear all
clc
figure
hold on
liminf=-12;
limsup=2;
liminf=-12;
limsup=2;


load ('CubeLinearResiduallambda1mu1Ceq1000Casym0Csmoothing5Dofs78C1F1.mat')
residual=residual(1:min(length(residual),50));
resSquare{1}=residual;
load ('CubeLinearResiduallambda1mu1Ceq1000Casym0Csmoothing5Dofs441C1F2.mat')
residual=residual(1:min(length(residual),50));
resSquare{2}=residual;
load ('CubeLinearResiduallambda1mu1Ceq1000Casym0Csmoothing5Dofs2967C1F3.mat')
residual=residual(1:min(length(residual),50));
resSquare{3}=residual;
load ('CubeLinearResiduallambda1mu1Ceq1000Casym0Csmoothing5Dofs21771C1F4.mat')
residual=residual(1:min(length(residual),50));
resSquare{4}=residual;
load ('CubeLinearResiduallambda1mu1Ceq1000Casym0Csmoothing5Dofs166803C1F5.mat')
residual=residual(1:min(length(residual),50));
resSquare{5}=residual;



plot(log10(resSquare{1}),'-b*','LineWidth',5,'MarkerSize',30);
plot(log10(resSquare{2}),'-g*','LineWidth',5,'MarkerSize',30);
plot(log10(resSquare{3}),'-r*','LineWidth',5,'MarkerSize',30);
plot(log10(resSquare{4}),'-k*','LineWidth',5,'MarkerSize',30);
plot(log10(resSquare{5}),'-c*','LineWidth',5,'MarkerSize',30);




l=legend('Dofs=78','Dofs=441','Dofs=2967','Dofs=21771','Dofs=166803');l.FontSize=46;
l.FontSize=60;
l.Location='northeast';
xx=xlabel('Iterations');
yy=ylabel('log10(Residual)');
set(gca,'fontsize',60)
title('C_{eq}=1e3, Smoothing=5')

ylim([liminf limsup])

figure
hold on
for ii=1:length(resSquare)
    for kk=1:length(resSquare{ii})-1
        if(log10(resSquare{ii}(kk+1))>liminf)
        rate{ii}(kk)=resSquare{ii}(kk+1)/resSquare{ii}(kk);
        end
    end
end

plot(rate{1},'-b*','LineWidth',5,'MarkerSize',30);
plot(rate{2},'-g*','LineWidth',5,'MarkerSize',30);
plot(rate{3},'-r*','LineWidth',5,'MarkerSize',30);
plot(rate{4},'-k*','LineWidth',5,'MarkerSize',30);
plot(rate{5},'-c*','LineWidth',5,'MarkerSize',30);

l2=legend('Dofs=78','Dofs=441','Dofs=2967','Dofs=21771','Dofs=166803');
l2.FontSize=46;

l2.Location='southeast';
xx=xlabel('Iterations');
yy=ylabel('Rate');
set(gca,'fontsize',60)
title('C_{eq}=1e3, Smoothing=5')

