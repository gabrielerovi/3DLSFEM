close all
clear all
clc
figure
hold on
liminf=-9.5;
limsup=2;
liminf=-9;
limsup=2;
% load ('SignoriniMINResiduallambda1e+50mu1Ceq100Casym0Csmoothing5.mat')
% residual=residual(1:min(length(residual),50));
% resSquare{3}=residual;
% load ('SignoriniMEANResiduallambda1e+50mu1Ceq100Casym0Csmoothing5.mat')
% residual=residual(1:min(length(residual),50));
% resSquare{2}=residual;
% load ('SignoriniResiduallambda1e+50mu1Ceq100Casym0Csmoothing5.mat')
% residual=residual(1:min(length(residual),50));
% resSquare{1}=residual;
% load ('SignoriniNoConstraintResiduallambda1e+50mu1Ceq100Casym0Csmoothing5.mat')
% residual=residual(1:min(length(residual),50));
% resSquare{4}=residual;
% 
% load ('SignoriniMINResiduallambda1mu1Ceq100Casym0Csmoothing5.mat')
% residual=residual(1:min(length(residual),50));
% resSquare{7}=residual;
% load ('SignoriniMEANResiduallambda1mu1Ceq100Casym0Csmoothing5.mat')
% residual=residual(1:min(length(residual),50));
% resSquare{6}=residual;
% load ('SignoriniResiduallambda1mu1Ceq100Casym0Csmoothing5.mat')
% residual=residual(1:min(length(residual),50));
% resSquare{5}=residual;
% load ('SignoriniNoConstraintResiduallambda1mu1Ceq100Casym0Csmoothing5.mat')
% residual=residual(1:min(length(residual),50));
% resSquare{8}=residual;



load ('SignoriniMINResiduallambda1e+50mu1Ceq100Casym0Csmoothing5.mat')
residual=residual(1:min(length(residual),50));
resSquare{3}=residual;
load ('SignoriniMEANResiduallambda1e+50mu1Ceq100Casym0Csmoothing5.mat')
residual=residual(1:min(length(residual),50));
resSquare{2}=residual;
load ('SignoriniResiduallambda1e+50mu1Ceq100Casym0Csmoothing5.mat')
residual=residual(1:min(length(residual),50));
resSquare{1}=residual;
load ('SignoriniMINResiduallambda1mu1Ceq100Casym0Csmoothing5.mat')
residual=residual(1:min(length(residual),50));
resSquare{6}=residual;
load ('SignoriniMEANResiduallambda1mu1Ceq100Casym0Csmoothing5.mat')
residual=residual(1:min(length(residual),50));
resSquare{5}=residual;
load ('SignoriniResiduallambda1mu1Ceq100Casym0Csmoothing5.mat')
residual=residual(1:min(length(residual),50));
resSquare{4}=residual;

plot(log10(resSquare{1}),'-b*','LineWidth',5,'MarkerSize',30);
plot(log10(resSquare{2}),'-g*','LineWidth',5,'MarkerSize',30);
plot(log10(resSquare{3}),'-r*','LineWidth',5,'MarkerSize',30);
plot(log10(resSquare{4}),'-bo','LineWidth',5,'MarkerSize',30);
plot(log10(resSquare{5}),'-go','LineWidth',5,'MarkerSize',30);
plot(log10(resSquare{6}),'-ro','LineWidth',5,'MarkerSize',30);

% plot(log10(resSquare{1}),'-b*','LineWidth',5,'MarkerSize',30);
% plot(log10(resSquare{2}),'-g*','LineWidth',5,'MarkerSize',30);
% plot(log10(resSquare{3}),'-r*','LineWidth',5,'MarkerSize',30);
% plot(log10(resSquare{4}),'-k*','LineWidth',5,'MarkerSize',30);
% plot(log10(resSquare{5}),'-bo','LineWidth',5,'MarkerSize',30);
% plot(log10(resSquare{6}),'-go','LineWidth',5,'MarkerSize',30);
% plot(log10(resSquare{7}),'-ro','LineWidth',5,'MarkerSize',30);
% plot(log10(resSquare{8}),'-ko','LineWidth',5,'MarkerSize',30);



l=legend('\lambda=\infty, \mu=1, a)','\lambda=\infty, \mu=1, b)','\lambda=\infty, \mu=1, c)',...
         '\lambda=1,  \mu=1, a)', '\lambda=1,  \mu=1, b)', '\lambda=1,  \mu=1, c)');
% l=legend('\lambda=\infty, \mu=1, a)','\lambda=\infty, \mu=1, b)','\lambda=\infty, \mu=1, c)','\lambda=\infty, \mu=1, d)',...
%          '\lambda=1,  \mu=1, a)', '\lambda=1,  \mu=1, b)', '\lambda=1,  \mu=1, c)','\lambda=\infty, \mu=1, d)');l.FontSize=46;
l.FontSize=60;
l.Location='northeast';
xx=xlabel('Iterations');
yy=ylabel('log10(Residual)');
set(gca,'fontsize',60)
title('C_{eq}=1e2, C_{compl}=1e1, Dofs=53186')

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
plot((rate{1}),'-b*','LineWidth',5,'MarkerSize',30);
plot((rate{2}),'-g*','LineWidth',5,'MarkerSize',30);
plot((rate{3}),'-r*','LineWidth',5,'MarkerSize',30);
plot((rate{4}),'-bo','LineWidth',5,'MarkerSize',30);
plot((rate{5}),'-go','LineWidth',5,'MarkerSize',30);
plot((rate{6}),'-ro','LineWidth',5,'MarkerSize',30);

% plot((rate{1}),'-b*','LineWidth',5,'MarkerSize',30);
% plot((rate{2}),'-g*','LineWidth',5,'MarkerSize',30);
% plot((rate{3}),'-r*','LineWidth',5,'MarkerSize',30);
% plot((rate{4}),'-k*','LineWidth',5,'MarkerSize',30);
% plot((rate{5}),'-bo','LineWidth',5,'MarkerSize',30);
% plot((rate{6}),'-go','LineWidth',5,'MarkerSize',30);
% plot((rate{7}),'-ro','LineWidth',5,'MarkerSize',30);
% plot((rate{8}),'-ko','LineWidth',5,'MarkerSize',30);

l=legend('\lambda=\infty, \mu=1, a)','\lambda=\infty, \mu=1, b)','\lambda=\infty, \mu=1, c)',...
         '\lambda=1,  \mu=1, a)', '\lambda=1,  \mu=1, b)', '\lambda=1,  \mu=1, c)');
% l=legend('\lambda=\infty, \mu=1, a)','\lambda=\infty, \mu=1, b)','\lambda=\infty, \mu=1, c)','\lambda=\infty, \mu=1, d)',...
%          '\lambda=1,  \mu=1, a)', '\lambda=1,  \mu=1, b)', '\lambda=1,  \mu=1, c)','\lambda=\infty, \mu=1, d)');l.FontSize=46;
l.FontSize=60;
l.Location='southeast';
xx=xlabel('Iterations');
yy=ylabel('Rate');
set(gca,'fontsize',60)
title('C_{eq}=1e2, C_{compl}=1e1, Dofs=53186')

