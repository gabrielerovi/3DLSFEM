function [x,norm_res] = ArnoldVcycleUTOPIA(A,b,x,vertices_per_processor,smoothing_steps,mesh_parameters,graph,Pnt,maxlev,lev)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Definition of the lev-Matrix with boundary conditions (also frictionless) Abc %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 F=lev;
 C=lev-1;
 RemoveF=mesh_parameters{F}.RemoveNT;
 Abc(RemoveF,:)=0;
 Abc(:, mesh_parameters{F}.E_contact+mesh_parameters{F}.NE)=0;
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % we assume b has only bc, not external forces <----------- not sure of
 % this, I am removing just bc dofs
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% On coarse levels, all the boundary conditions are homogeneous %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if(maxlev>lev)
 Abc(:,RemoveF)=0;
 b(RemoveF)=0;
 end
 for ii=RemoveF
 Abc(ii,ii)=1;
 end

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% If on coarsest level: compute activeset %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(F==1)
    x=Abc\b;
else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Otherwise: 1) pre-smoothing, 2) project for new correction and 3) post-smoothing %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Lev>1: pre-smoothing %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    x=UTOPIASmoother(A,x,b,vertices_per_processor{F},mesh_parameters{F}.Patch_Internal_All,smoothing_steps); 
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Fine residual and coarse problem (AC correctionC = resC, resC <= ConstraintC) %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    resF= b - Abc * x;  
    resF(RemoveF)=0;
    resC=P'*resF;
    resC(RemoveC)=0;
    corrC=sparse(length(resC),1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Compute the correction correctionC: AC correctionC = resC, resC <= ConstraintC %%%%
    %%%%% Then compute the fine correction and the actual solution x                     %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
     corrC = ArnoldVcycleUTOPIA(ACnobc,resC,corrC,vertices_per_processor,smoothing_steps,Patch_Internal_All,graph,Pnt,maxlev,C)
     corrF=P*corrC;
     corrF(RemoveF)=0;    
     x=x+corrF;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Lev>1: pre-smoothing %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x=UTOPIASmoother(A,x,b,vertices_per_processor{F},mesh_parameters{F}.Patch_Internal_All,smoothing_steps); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Computation of the norm of the residual, without considering bcs and active set dofs %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(lev==maxlev)
    resF=b -Abc*x;
    resF(RemoveF)=0;
    norm_res=norm(resF);
end

end

