function [x,WorkingSet,norm_res,norm_corrF] = ArnoldVcycleContact3(Complementarity,h,A,b,x,Constraint,WorkingSet,parameters,mesh,maps,graph,Pnt,maxlev,lev)
 % FINE && COARSE %
 F=lev;
 C=lev-1;
 NEF=mesh{F}.NE;
 NF=mesh{F}.N;
 RemoveF=mesh{F}.RemoveNT;
 

 % A is without bc, Abc is with (also with frictionless)
 
 Aenergy=A+ h(lev)/h(maxlev)*Complementarity;
 Abc=Aenergy;
 Abc(RemoveF,:)=0;
 Abc(:, mesh{F}.E_contact+mesh{F}.NE)=0;
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % we assume b has only bc, not external forces
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % on coarse levels, all the boundary conditions are homogeneous 
 if(maxlev>lev)
 Abc(:,RemoveF)=0;
 b(RemoveF)=0;
 end

 for ii=RemoveF
 Abc(ii,ii)=1;
 end

 P=Pnt{ max(1,C)};

if(F==1)
    norm_corrF=0;
    % COARSE LEVEL: ACTIVE SET
    norm_res=[];
    WorkingSet_Loc=[];
    B=speye(length(b));
    [x,lambda,WorkingSet_Loc] = ArnoldActiveset2(Abc,B,b,Constraint,WorkingSet_Loc);
%      [x,WorkingSet]=activesetResidual2(Abc,b,Constraint,WorkingSet_Loc);
     
else
    % LEV>1: SMOOTHING
    NEC=mesh{C}.NE;
    NC=mesh{C}.N;
    RemoveC=mesh{C}.RemoveNT;
    smoothing_steps=parameters.smoothing_steps;
    is_on_coarser_grid=0;
    top2bottom=1;

    % PRE-SMOOTHING
    [x,WorkingSet,energy]=ArnoldSmootherContact7(WorkingSet,Aenergy,Abc,x,b,Constraint,mesh{F},maps.Patch_Internal_All{F},smoothing_steps,graph{F},top2bottom) ;
    ActiveSet=find(WorkingSet>0);

    % TRUNCATED BASES
    P(ActiveSet,:)=0;

    % COARSE GAP FUNCTION
    ConstraintF=Constraint-x;
    ConstraintF(ActiveSet)=10^10;
    ConstraintC= ArnoldCoarseConstraint(ConstraintF,mesh,maps,C,F);

    % FINE RESIDUAL
    resF= b - Abc * x;  
    resF(RemoveF)=0;

    % COARSE TRUNCATED RESIDUAL
    resC=P'*resF;
    resC(RemoveC)=0;

    % COARSE TRUNCATED MATRIX (homogeneous bc)
    ACnobc=P'*A*P;
    ComplementarityC=P'* Complementarity  *P;
    AC=ACnobc;
    AC(RemoveC,:)=0;
    AC(:,RemoveC)=0;
    for ii=RemoveC
        AC(ii,ii)=1;
    end

    
    corrC=sparse(2*NEC+2*NC,1);
    WorkingSetC=sparse(length(corrC),1);
     [corrC,WorkingSetC,norm_res,norm_corrF] = ArnoldVcycleContact3(ComplementarityC,h,ACnobc,resC,corrC,ConstraintC,WorkingSetC,parameters,mesh,maps,graph,Pnt,maxlev,C);

    corrF=P*corrC;
    corrF(RemoveF)=0;
    
    v1=0.5*x'*Aenergy*x-b'*x
    x=x+corrF;
    v2=0.5*x'*Aenergy*x-b'*x
    v2<v1
    constraint_ok=sum(x-Constraint>10^(-12));
    if(constraint_ok>0)
        fermami=1
    end
    norm_corrF=norm(corrF);
    % POST-SMOOTHING
    top2bottom=0;
    [x,WorkingSet,energy]=ArnoldSmootherContact7(WorkingSet,Aenergy,Abc,x,b,Constraint,mesh{F},maps.Patch_Internal_All{F},smoothing_steps,graph{F},top2bottom) ;

resF=b -Abc*x;
resF(RemoveF)=0;
resF(find(WorkingSet==1))=0;

norm_res=norm(resF);


end


end






% [M_Normal_Tangent,M_Normal_TangentT] = MatrixOnGammaCwithNormalTangentComponents(mesh{C},1);
% sol=M_Normal_Tangent*corrC;
% print_displacement_solution2(C,mesh,sol(1+2*mesh{C}.NE:2*mesh{C}.NE+mesh{C}.N)',sol(1+2*mesh{C}.NE+mesh{C}.N:end)');
% 
% [M_Normal_Tangent,M_Normal_TangentT] = MatrixOnGammaCwithNormalTangentComponents(mesh{F},1);
% sol=M_Normal_Tangent*corrF;
% print_displacement_solution2(F,mesh,sol(1+2*mesh{F}.NE:2*mesh{F}.NE+mesh{F}.N)',sol(1+2*mesh{F}.NE+mesh{F}.N:end)');
% 
% [M_Normal_Tangent,M_Normal_TangentT] = MatrixOnGammaCwithNormalTangentComponents(mesh{F},1);
% sol=M_Normal_Tangent*x;
% print_displacement_solution2(F,mesh,sol(1+2*mesh{F}.NE:2*mesh{F}.NE+mesh{F}.N)',sol(1+2*mesh{F}.NE+mesh{F}.N:end)');

