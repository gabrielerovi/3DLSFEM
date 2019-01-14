function [ERRUX ,ERRUY, ERRSXX, ERRSXY, ERRSYX, ERRSYY, ERRSXXdisp,ERRSXYdisp, ERRSYXdisp, ERRSYYdisp]=LS_error_computation(parameters,mesh,x)
   % L2ErrDisp,H1ErrDisp,L2ErrStress,HdivErrStress]=LS_error_computation(mesh,x,DispX,DispY,StressXX,StressXY,StressYX,StressYY)

mu=parameters.mu;
lambda=parameters.lambda;
DispX=parameters.disp1;
DispY=parameters.disp2;
StressXX=parameters.stressxx;
StressXY=parameters.stressxy;
StressYX=parameters.stressyx;
StressYY=parameters.stressyy;

L=size(mesh);
L=L(1);

N=mesh{L}.N;
NE=mesh{L}.NE;
NT=mesh{L}.NT;
 if(strcmp(parameters.input_name  ,'DispElasticity'))
sigmax=zeros(NE,1);
sigmay=zeros(NE,1);
ux=    x(1:N);
uy=    x(1+N:end);
 else
sigmax=x(1:NE);
sigmay=x(1+NE:2*NE);
ux=    x(1+2*NE:2*NE+N);
uy=    x(1+2*NE+N:end);
 end


qrule=5;
    L2ErrDisp    = 0;
    H1ErrDisp    = 0;
    L2ErrStress  = 0;
    HdivErrStress = 0;
    ERRUX  = 0;
    ERRUY  = 0;
    
    ERRSXX = 0;
    ERRSXY = 0;
    ERRSYX = 0;
    ERRSYY = 0;
    
    ERRSXXdisp = 0;
    ERRSXYdisp = 0;
    ERRSYXdisp = 0;
    ERRSYYdisp = 0;
    

    for t=1:mesh{L}.NT
    
    % nodal and edge dofs
    elem=mesh{L}.elem(t,:);
    elemE=mesh{L}.elemE(t,:);
    node=mesh{L}.node(elem,:);
    
    % quadrature points
    [q_point,weights,area]=quadrature_points_2D(qrule,node);
    number_of_qp=length(q_point)
    % local ux,uy,sigmax,sigmay
    uxLoc=ux(elem);
    uyLoc=uy(elem);
    sigmaxLoc=sigmax(elemE);
    sigmayLoc=sigmay(elemE);  
    
    % P1 and RT0 shape functions
    [RT0,RT0_div] = phiRT2Dcell(q_point,node);
     P1_basis = phiP12D(q_point,node);
     P1grad = P1grad2D(node);
     
         L2ErrDispLoc = 0;
         H1ErrDispLoc = 0;
         L2ErrStressLoc = 0;
         HdivErrStressLoc = 0;
         
         for qp=1:number_of_qp
         
             UXLOC=0;
             UYLOC=0;
             SXXLOC=0;
             SXYLOC=0;
             SYXLOC=0;
             SYYLOC=0;
             
             sxxloc=0;
             sxyloc=0;
             syxloc=0;
             syyloc=0;
             
             for nn=1:3
                 UXLOC = UXLOC + uxLoc(nn)*P1_basis(nn,qp) ;
                 UYLOC = UYLOC + uyLoc(nn)*P1_basis(nn,qp) ;
                 SXXLOC = SXXLOC + sigmaxLoc(nn) * RT0{nn}(qp,1);
                 SXYLOC = SXYLOC + sigmaxLoc(nn) * RT0{nn}(qp,2);
                 SYXLOC = SYXLOC + sigmayLoc(nn) * RT0{nn}(qp,1);
                 SYYLOC = SYYLOC + sigmayLoc(nn) * RT0{nn}(qp,2);  
                 
                 sxxloc= sxxloc + 2 * mu *(uxLoc(nn)*P1grad(nn,1) ) + lambda * (uxLoc(nn)*P1grad(nn,1) + uyLoc(nn)*P1grad(nn,2) );
                 sxyloc= sxyloc + 2 * mu * 0.5 * (uxLoc(nn)*P1grad(nn,2) + uyLoc(nn)*P1grad(nn,1));
                 syxloc= syxloc + 2 * mu * 0.5 * (uxLoc(nn)*P1grad(nn,2) + uyLoc(nn)*P1grad(nn,1));
                 syyloc= syyloc + 2 * mu * (uyLoc(nn)*P1grad(nn,2) )+ lambda * (uxLoc(nn)*P1grad(nn,1) + uyLoc(nn)*P1grad(nn,2) );
             end
             qx=q_point(qp,1);
             qy=q_point(qp,2);
           ERRUX  = ERRUX + area * weights(qp) * (UXLOC - DispX(qx,qy))^2;
           ERRUY  = ERRUY + area * weights(qp) * (UYLOC - DispY(qx,qy))^2;
           ERRSXX = ERRSXX+ area * weights(qp) * (SXXLOC - StressXX(qx,qy))^2;
           ERRSXY = ERRSXY+ area * weights(qp) * (SXYLOC - StressXY(qx,qy))^2;
           ERRSYX = ERRSYX+ area * weights(qp) * (SYXLOC - StressYX(qx,qy))^2;
           ERRSYY = ERRSYY+ area * weights(qp) * (SYYLOC - StressYY(qx,qy))^2;
           
           ERRSXXdisp = ERRSXXdisp + area * weights(qp) * (sxxloc - StressXX(qx,qy))^2;
           ERRSXYdisp = ERRSXYdisp + area * weights(qp) * (sxyloc - StressXY(qx,qy))^2;
           ERRSYXdisp = ERRSYXdisp + area * weights(qp) * (syxloc - StressYX(qx,qy))^2;
           ERRSYYdisp = ERRSYYdisp + area * weights(qp) * (syyloc - StressYY(qx,qy))^2;
         end
     

    end
           ERRUX  = sqrt(ERRUX);
           ERRUY  = sqrt(ERRUY);
           
           ERRSXX = sqrt(ERRSXX);
           ERRSXY = sqrt(ERRSXY);
           ERRSYX = sqrt(ERRSYX);
           ERRSYY = sqrt(ERRSYY);
           
           ERRSXXdisp = sqrt(ERRSXXdisp);
           ERRSXYdisp = sqrt(ERRSXYdisp);
           ERRSYXdisp = sqrt(ERRSYXdisp);
           ERRSYYdisp = sqrt(ERRSYYdisp);
           
 
end






