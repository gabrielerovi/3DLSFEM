function [erruL2,errgrad,errsigmaL2,errdivsigmaL2]=error_displacement(mesh,x,parameters)
L=size(mesh);
L=L(1);
figure
increase=1;
cont_elem=1;
max_err=0;
node_per_elem=mesh{L}.node_per_elem;

disp1ex=parameters.disp1;
disp2ex=parameters.disp2;
disp3ex=parameters.disp3;

dispgradx=parameters.dispgradx;
dispgrady=parameters.dispgrady;
dispgradz=parameters.dispgradz;

stressxx=parameters.stressxx; stressxy=parameters.stressxy; stressxz=parameters.stressxz;
stressyx=parameters.stressyx; stressyy=parameters.stressyy; stressyz=parameters.stressyz;
stresszx=parameters.stresszx; stresszy=parameters.stresszy; stresszz=parameters.stresszz;

divstress1=parameters.divstress1;
divstress2=parameters.divstress2;
divstress3=parameters.divstress3;

NF=mesh{L}.NF;
NF3=mesh{L}.NF*3;
N=mesh{L}.N;

sigmatot1=x(1:NF);
sigmatot2=x(1+NF:2*NF);
sigmatot3=x(1+2*NF:NF3);

disp1=x(1+NF3:NF3+N);
disp2=x(1+N+NF3:NF3+2*N);
disp3=x(1+2*N+NF3:NF3+3*N);

hold on

err1=zeros(N,1);
err2=zeros(N,1);
err3=zeros(N,1);
noderef=[0     0     0;
      1     0     0;
      0     1     0;
      0     0     1];
  qrule=4;
[q_pointref,weights,Volume]=quadrature_points_3D(qrule,noderef);
maxerr1=0;maxerr2=0;maxerr3=0;
for nn=1:mesh{L}.N
    node=mesh{L}.node(nn,1:parameters.dim);
    
    err1(nn,1)=(disp1(nn)-disp1ex(node(1),node(2),node(3) ))^2;
    err2(nn,1)=(disp1(nn)-disp2ex(node(1),node(2),node(3) ))^2;
    err3(nn,1)=(disp1(nn)-disp3ex(node(1),node(2),node(3) ))^2;
    
    maxerr1=max(maxerr1,abs(disp1(nn)-disp1ex(node(1),node(2),node(3) )));
    maxerr2=max(maxerr2,abs(disp2(nn)-disp2ex(node(1),node(2),node(3) )));
    maxerr3=max(maxerr3,abs(disp3(nn)-disp3ex(node(1),node(2),node(3) )));
end

err1=0; errgrad1=0; errsigma1=0; errdivsigma1=0;
err2=0; errgrad2=0; errsigma2=0; errdivsigma2=0;
err3=0; errgrad3=0; errsigma3=0; errdivsigma3=0;

for tt=1:mesh{L}.NT
    elem=mesh{L}.elem(tt,:);
    elemF=mesh{L}.elemF(tt,:);
    node=mesh{L}.node(elem,1:3);
    [J,detJ]=J_and_detJ(node); 
    q_point=q_pointref*J';
    absdetJ=abs(detJ);
    [P1_basis,P1grad] = phiP13D(q_point,node);
    disp1exqp=disp1ex(q_point(:,1),q_point(:,2),q_point(:,3));
    disp2exqp=disp2ex(q_point(:,1),q_point(:,2),q_point(:,3));
    disp3exqp=disp3ex(q_point(:,1),q_point(:,2),q_point(:,3));
    
    grad1exqp=dispgradx(q_point(:,1),q_point(:,2),q_point(:,3));
    grad2exqp=dispgrady(q_point(:,1),q_point(:,2),q_point(:,3));
    grad3exqp=dispgradz(q_point(:,1),q_point(:,2),q_point(:,3));

    disp1loc=disp1(elem);
    disp2loc=disp2(elem);
    disp3loc=disp3(elem);

    grad1=sum(disp1loc.*P1grad);grad1=repmat(grad1,length(q_point(:,1)),1); gradxx=grad1(:,1);gradxy=grad1(:,2);gradxz=grad1(:,3);
    grad2=sum(disp2loc.*P1grad);grad2=repmat(grad2,length(q_point(:,1)),1); gradyx=grad2(:,1);gradyy=grad2(:,2);gradyz=grad2(:,3);
    grad3=sum(disp3loc.*P1grad);grad3=repmat(grad3,length(q_point(:,1)),1); gradzx=grad3(:,1);gradzy=grad3(:,2);gradzz=grad3(:,3);
    
    
    
    
    disp1qp=sum(disp1loc'.*P1_basis',2);
    disp2qp=sum(disp2loc'.*P1_basis',2);
    disp3qp=sum(disp3loc'.*P1_basis',2);
    
    err1=err1+absdetJ *weights.*(disp1qp-disp1exqp).^2;
    err2=err2+absdetJ *weights.*(disp1qp-disp1exqp).^2;
    err3=err3+absdetJ *weights.*(disp1qp-disp1exqp).^2;
    
    errgrad1=errgrad1+absdetJ *weights.*( (gradxx-grad1exqp).^2 + (gradxy-grad2exqp).^2 + (gradxz-grad3exqp).^2);
    errgrad2=errgrad2+absdetJ *weights.*( (gradyx-grad1exqp).^2 + (gradyy-grad2exqp).^2 + (gradyz-grad3exqp).^2);
    errgrad3=errgrad3+absdetJ *weights.*( (gradzx-grad1exqp).^2 + (gradzy-grad2exqp).^2 + (gradzz-grad3exqp).^2);
    
    
    sigma1loc=sigmatot1(elemF);
    sigma2loc=sigmatot2(elemF);
    sigma3loc=sigmatot3(elemF);
    
    [RT0_basis,RT0_divergence] = phiRT3Dcell(q_point,node);
    
    sigma1=sigma1loc(1)*RT0_basis{1}+sigma1loc(2)*RT0_basis{2}+sigma1loc(3)*RT0_basis{3}+sigma1loc(4)*RT0_basis{4};
    sigma2=sigma2loc(1)*RT0_basis{1}+sigma2loc(2)*RT0_basis{2}+sigma2loc(3)*RT0_basis{3}+sigma2loc(4)*RT0_basis{4};
    sigma3=sigma3loc(1)*RT0_basis{1}+sigma3loc(2)*RT0_basis{2}+sigma3loc(3)*RT0_basis{3}+sigma3loc(4)*RT0_basis{4};
    
    sigmaxx=sigma1(:,1); sigmaxy=sigma1(:,2); sigmaxz=sigma1(:,3);
    sigmayx=sigma2(:,1); sigmayy=sigma2(:,2); sigmayz=sigma2(:,3);
    sigmazx=sigma3(:,1); sigmazy=sigma3(:,2); sigmazz=sigma3(:,3);
    
    sigmaxxex=stressxx(q_point(:,1),q_point(:,2),q_point(:,3)); sigmaxyex=stressxy(q_point(:,1),q_point(:,2),q_point(:,3)); sigmaxzex=stressxz(q_point(:,1),q_point(:,2),q_point(:,3));
    sigmayxex=stressyx(q_point(:,1),q_point(:,2),q_point(:,3)); sigmayyex=stressyy(q_point(:,1),q_point(:,2),q_point(:,3)); sigmayzex=stressyz(q_point(:,1),q_point(:,2),q_point(:,3));
    sigmazxex=stresszx(q_point(:,1),q_point(:,2),q_point(:,3)); sigmazyex=stresszy(q_point(:,1),q_point(:,2),q_point(:,3)); sigmazzex=stresszz(q_point(:,1),q_point(:,2),q_point(:,3));

    
    errsigma1=errsigma1+absdetJ *weights.*( (sigmaxx-sigmaxxex).^2 + (sigmaxy-sigmaxyex).^2 + (sigmaxz-sigmaxzex).^2);
    errsigma2=errsigma2+absdetJ *weights.*( (sigmayx-sigmayxex).^2 + (sigmayy-sigmayyex).^2 + (sigmayz-sigmayzex).^2);
    errsigma3=errsigma3+absdetJ *weights.*( (sigmazx-sigmazxex).^2 + (sigmazy-sigmazyex).^2 + (sigmazz-sigmazzex).^2);

    
    divsigma1ex=divstress1(q_point(:,1),q_point(:,2),q_point(:,3));
    divsigma2ex=divstress2(q_point(:,1),q_point(:,2),q_point(:,3));
    divsigma3ex=divstress3(q_point(:,1),q_point(:,2),q_point(:,3));
    
    
    RT0_divergence1=RT0_divergence*sigma1loc; RT0_divergence1=repmat(RT0_divergence1,length(q_point(:,1)),1);
    RT0_divergence2=RT0_divergence*sigma2loc; RT0_divergence2=repmat(RT0_divergence2,length(q_point(:,1)),1);
    RT0_divergence3=RT0_divergence*sigma3loc; RT0_divergence3=repmat(RT0_divergence3,length(q_point(:,1)),1);
    
    errdivsigma1=errdivsigma1+absdetJ *weights.*(RT0_divergence1-divsigma1ex).^2 ;
    errdivsigma2=errdivsigma2+absdetJ *weights.*(RT0_divergence2-divsigma2ex).^2 ;
    errdivsigma3=errdivsigma3+absdetJ *weights.*(RT0_divergence3-divsigma3ex).^2 ;

    
end

erruL2=sqrt(sum(err1)+sum(err2)+sum(err3));
errgrad=sqrt(sum(errgrad1)+sum(errgrad2)+sum(errgrad3));
errsigmaL2=sqrt(sum(errsigma1)+sum(errsigma2)+sum(errsigma3));
errdivsigmaL2=sqrt(sum(errdivsigma1)+sum(errdivsigma2)+sum(errdivsigma3));
end