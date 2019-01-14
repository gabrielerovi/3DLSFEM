function [x,WorkingSet,norm_res] = ArnoldVcycleContact3DEdge(A,Complementarity,b,x,h,Constraint,WorkingSet,smoothing_steps,mesh_parameters,Pnt,maxlev,lev)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Definition of the lev-Matrix with boundary conditions (also frictionless) Abc %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 F=lev;
 C=lev-1;
 RemoveF=mesh_parameters{F}.RemoveNT;
 Abc=A+ (h(lev)/h(maxlev))^2*Complementarity;
%  Abc(RemoveF,:)=0;
%  Abc(:, mesh_parameters{F}.E_contact+mesh_parameters{F}.NF)=0;
 
 % RICONTROLLAAAA........................................................................................................................................................
  % RICONTROLLAAAA........................................................................................................................................................
  
  %%% RemoveNT???????????????????????????????????????????????????????????????????????????
%  Abc=add_boundary_bc_system(Abc,mesh_parameters{F}.RemoveNT);
 [Abc,b]=add_boundary_bc_system_symmetric(A,b,mesh_parameters{F}.RemoveNT);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % we assume b has only bc, not external forces <----------- not sure of
 % this, I am removing just bc dofs
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% On coarse levels, all the boundary conditions are homogeNFous %%%%%
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
    WorkingSet_Loc=[];
    B=speye(length(b));
    [x,lambda,WorkingSet] = ArnoldActiveset(Abc,B,b,Constraint,WorkingSet_Loc);   
   % [x,lambda,WorkingSet] =  ActiveSet(Abc,B,b,Constraint,WorkingSet_Loc,x)  ;   
else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Otherwise: 1) pre-smoothing, 2) projFCt for NFw corrFCtion and 3) post-smoothing %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Lev>1: pre-smoothing %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    [x,WorkingSet]=ArnoldSmootherContactEdge(WorkingSet,Abc,x,b,Constraint,flip_vector(mesh_parameters{F}.edgegraphinterior,0),flip_vector(mesh_parameters{F}.edgegraphgammac,0), mesh_parameters{F}.Edge_Patch_Internal_All,smoothing_steps) ;
    ActualActiveSet=find(WorkingSet>0);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Coarse gap function %%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ConstraintF=Constraint-x;
    ConstraintF(ActualActiveSet)=10^10;
    
    mp.NFC=mesh_parameters{C}.NF;
    mp.NFF=mesh_parameters{F}.NF;
    mp.NC=mesh_parameters{C}.N;
    mp.NF=mesh_parameters{F}.N;
    mp.FC_contact=mesh_parameters{C}.F_contact;
    mp.Patch_Node_Monotone=mesh_parameters{C}.Patch_Node_Monotone;
    mp.Patch_Face_Monotone=mesh_parameters{C}.Patch_Face_Monotone;
    ConstraintC= ArnoldCoarseConstraint3D(ConstraintF,mp);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Truncated basis: remove all the contributions of dofs on the ActiveSet %%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    P=Pnt{C};
    P(ActualActiveSet,:)=0;
    RemovFC=mesh_parameters{C}.RemoveNT;

     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% FiNF residual and coarse problem (AC corrFCtionC = resC, resC <= ConstraintC) %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    resF= b - Abc * x;  
    resF(RemoveF)=0;
    resC=P'*resF;
    resC(RemovFC)=0;
    ACnobc=P'*A*P;
    ComplementarityC=P'* Complementarity  *P;
    corrC=sparse(length(resC),1);
    WorkingSetC=sparse(length(resC),1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Compute the corrFCtion corrFCtionC: AC corrFCtionC = resC, resC <= ConstraintC %%%%
    %%%%% Then compute the fiNF corrFCtion and the actual solution x                     %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    [corrC,WorkingSetC] = ArnoldVcycleContact3DEdge(ACnobc,ComplementarityC,resC,corrC,h,ConstraintC,WorkingSetC,smoothing_steps,mesh_parameters,Pnt,maxlev,C);
     corrF=P*corrC;
     corrF(RemoveF)=0;    
      x=x+corrF;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Lev>1: pre-smoothing %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [x,WorkingSet]=ArnoldSmootherContactEdge(WorkingSet,Abc,x,b,Constraint,flip_vector(mesh_parameters{F}.edgegraphinterior,1),flip_vector(mesh_parameters{F}.edgegraphgammac,1), mesh_parameters{F}.Edge_Patch_Internal_All,smoothing_steps) ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Computation of the norm of the residual, without considering bcs and active set dofs %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(lev==maxlev && maxlev>1)
    resF=b -Abc*x;
    resF(RemoveF)=0;
    resF(find(WorkingSet==1))=0;
    norm_res=norm(resF);
elseif(lev==maxlev && maxlev==1)
     resF=b -Abc*x;
    resF(RemoveF)=0;
    resF(WorkingSet)=0;
    norm_res=norm(resF);   
end

end
