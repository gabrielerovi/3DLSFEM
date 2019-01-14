close all
clear all
figure
hold on

liminf=-9.5;
limsup=-2;


load ('CubeResiduallambda1mu1Ceq1000Casym0Csmoothing3Dofs441C1F2.mat')
resSquare{1}=residual;
load ('CubeResiduallambda1mu1Ceq1000Casym0Csmoothing3Dofs2967C1F3.mat')
resSquare{2}=residual;
load ('CubeResiduallambda1mu1Ceq1000Casym0Csmoothing3Dofs21771C1F4.mat')
resSquare{3}=residual;
load ('CubeResiduallambda1mu1Ceq1000Casym0Csmoothing3Dofs166803C1F5.mat')
resSquare{4}=residual;
load ('CubeResiduallambda1e+50mu1Ceq1000Casym0Csmoothing3Dofs441C1F2.mat')
resSquare{5}=residual;
load ('CubeResiduallambda1e+50mu1Ceq1000Casym0Csmoothing3Dofs2967C1F3.mat')
resSquare{6}=residual;
load ('CubeResiduallambda1e+50mu1Ceq1000Casym0Csmoothing3Dofs21771C1F4.mat')
resSquare{7}=residual;
load ('CubeResiduallambda1e+50mu1Ceq1000Casym0Csmoothing3Dofs166803C1F5.mat')
resSquare{8}=residual;

% for lev=1:length(resSquare)
% tmp=find(resSquare{lev}<-9);
% minimum(lev)=tmp(1);
% % colors(lev,1:3)=0.0001+lev/(length(resSquare)+1);
% % colors(lev,2)=0.0001+(length(resSquare)-lev)/(length(resSquare)+1);;
% % colors(lev,3)=0.0001+lev/(length(resSquare)+1);
% end
hold on
colors = {'k','b','r','g','k','b','r','g'};
for lev=1:length(resSquare)/2
plot(log10(resSquare{lev}),'color',colors{lev},'LineWidth',12);
end
for lev=1+length(resSquare)/2:length(resSquare)
plot(log10(resSquare{lev}),'color',colors{lev},'LineWidth',12,'Marker','o','MarkerSize',24,'MarkerFaceColor',[1 1 1]);
end

l=legend('L=2, Dofs=441, \lambda=1 ','L=3, Dofs=2967, \lambda=1 ','L=4, Dofs=21771, \lambda=1 ','L=5, Dofs=166803, \lambda=1 ', ...
         'L=2, Dofs=441, \lambda=\infty ','L=3, Dofs=2967, \lambda=\infty ','L=4, Dofs=21771, \lambda=\infty ','L=5, Dofs=166803, \lambda=\infty ');
l.FontSize=46;
xx=xlabel('Iterations');
yy=ylabel('log10(Residual)');
set(gca,'fontsize',48)
title('\lambda=1, \mu=1, C_{eq}=1e2, C_{compl}=1e2')

ylim([liminf limsup])
for ii=1:length(resSquare)
    for kk=1:length(resSquare{ii})   
        if(log10(resSquare{ii}(kk))<liminf+0.05)
        resSquare{ii}(kk:end)=[];
        break
        end
    end
end
    
for ii=1:8
    for kk=1:length(resSquare{ii})-1
%         if(kk==8 && ii==3)
%             fermami=1
%         end
        if((resSquare{ii}(kk))>liminf+0.5)
        rate{ii}(kk)=(resSquare{ii}(kk+1))/(resSquare{ii}(kk));
        end
    end
end

figure 
hold on

for lev=1:length(resSquare)/2
plot(rate{lev},'color',colors{lev},'LineWidth',12);
end
for lev=1+length(resSquare)/2:length(resSquare)
plot(rate{lev},'color',colors{lev},'LineWidth',12,'Marker','o','MarkerSize',24,'MarkerFaceColor',[1 1 1]);
end



l=legend('L=2, Dofs=441, \lambda=1 ','L=3, Dofs=2967, \lambda=1 ','L=4, Dofs=21771, \lambda=1 ','L=5, Dofs=166803, \lambda=1 ', ...
         'L=2, Dofs=441, \lambda=\infty ','L=3, Dofs=2967, \lambda=\infty ','L=4, Dofs=21771, \lambda=\infty ','L=5, Dofs=166803, \lambda=\infty ');
l.FontSize=46;
l.Location='southeast';
xx=xlabel('Iterations');
yy=ylabel('Rate');
set(gca,'fontsize',48)
title('\lambda=1, \mu=1, C_{eq}=1e2, C_{compl}=1e2')
