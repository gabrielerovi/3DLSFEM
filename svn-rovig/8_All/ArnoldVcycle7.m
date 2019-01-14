function [x] = ArnoldVcycle7(A,b,x,parameters,mesh,maps,graph,Pnt,maxlev,lev)
 % FINE && COARSE %
 F=lev;
 C=lev-1;
 NEF=mesh{F}.NE;
 NF=mesh{F}.N;
 RemoveF=mesh{F}.RemoveNT;
 

 % A is without bc, Abc is with (also with frictionless)
 Abc=A;
 Abc(RemoveF,:)=0;
 Abc(:, mesh{F}.E_contact+mesh{F}.NE)=0;
 % on coarse levels, all the boundary conditions are homogeneous 
 if(maxlev>lev)
 Abc(:,RemoveF)=0;
 bnt(RemoveF)=0;
 end

 for ii=RemoveF
 Abc(ii,ii)=1;
 end

 P=Pnt{ max(1,C)};

if(F==1)
    % COARSE LEVEL: ACTIVE SET
    norm_res=[];
    x=Abc\b;
else
    % LEV>1: SMOOTHING
    NEC=mesh{C}.NE;
    NC=mesh{C}.N;
    RemoveC=mesh{C}.RemoveNT;
    smoothing_steps=parameters.smoothing_steps;
    is_on_coarser_grid=0;
    top2bottom=1;

    % PRE-SMOOTHING
    x=ArnoldSmoother7(Abc,x,b,maps.Patch_Internal_All{F},smoothing_steps,graph{F},top2bottom) ;

    % FINE RESIDUAL
    resF= b - Abc * x;  
    resF(RemoveF)=0;

    % COARSE RESIDUAL
    resC=P'*resF;
    resC(RemoveC)=0;

    % COARSE MATRIX (homogeneous bc)
    ACnobc=P'*A*P;
    AC=ACnobc;
    AC(RemoveC,:)=0;
    AC(:,RemoveC)=0;
    for ii=RemoveC
        AC(ii,ii)=1;
    end

    
      corrC=sparse(2*NEC+2*NC,1);
      corrC = ArnoldVcycle7(ACnobc,resC,corrC,parameters,mesh,maps,graph,Pnt,maxlev,C);

      corrF=P*corrC;
      corrF(RemoveF)=0;
      x=x+corrF;
 
    % POST-SMOOTHING
    x=ArnoldSmoother7(Abc,x,b,maps.Patch_Internal_All{F},smoothing_steps,graph{F},top2bottom) ;
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


