function [x,WorkingSet,norm_res] = ArnoldVcycleContactNoCoarseConstraint(A,Complementarity,b,x,h,Constraint,WorkingSet,smoothing_steps,mesh_parameters,graph,Pnt,maxlev,lev)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Definition of the lev-Matrix with boundary conditions (also frictionless) Abc %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 F=lev;
 C=lev-1;
 RemoveF=mesh_parameters{F}.RemoveNT;
 Abc=A+ h(lev)/h(maxlev)*Complementarity;
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
    WorkingSet_Loc=[];
    B=speye(length(b));
    [x,lambda,WorkingSet_Loc] = ArnoldActiveset(Abc,B,b,Constraint,WorkingSet_Loc);     
else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Otherwise: 1) pre-smoothing, 2) project for new correction and 3) post-smoothing %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Lev>1: pre-smoothing %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    [x,WorkingSet]=ArnoldSmootherContact(WorkingSet,Abc,x,b,Constraint,flip_vector(graph{F},0),mesh_parameters{F}.N_contact,mesh_parameters{F}.Patch_Internal_All,smoothing_steps) ;
    ActiveSet=find(WorkingSet>0);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Truncated basis: remove all the contributions of dofs on the ActiveSet %%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    P=Pnt{C};
    P(ActiveSet,:)=0;
    RemoveC=mesh_parameters{C}.RemoveNT;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Coarse gap function %%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ConstraintF=Constraint-x;
    ConstraintF(ActiveSet)=10^10;
    
    mp.NEC=mesh_parameters{C}.NE;
    mp.NEF=mesh_parameters{F}.NE;
    mp.NC=mesh_parameters{C}.N;
    mp.NF=mesh_parameters{F}.N;
    mp.EC_contact=mesh_parameters{C}.E_contact;
    mp.Patch_Node_Monotone=mesh_parameters{C}.Patch_Node_Monotone;
    mp.Patch_Edge_Monotone=mesh_parameters{C}.Patch_Edge_Monotone;
    
    ConstraintC= 10^10 *ones(length(P(1,:)),1);


     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Fine residual and coarse problem (AC correctionC = resC, resC <= ConstraintC) %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    resF= b - Abc * x;  
    resF(RemoveF)=0;
    resC=P'*resF;
    resC(RemoveC)=0;
    ACnobc=P'*A*P;
    ComplementarityC=P'* Complementarity  *P;
    corrC=sparse(length(resC),1);
    WorkingSetC=sparse(length(resC),1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Compute the correction correctionC: AC correctionC = resC, resC <= ConstraintC %%%%
    %%%%% Then compute the fine correction and the actual solution x                     %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    [corrC,WorkingSetC] = ArnoldVcycleContactNoCoarseConstraint(ACnobc,ComplementarityC,resC,corrC,h,ConstraintC,WorkingSetC,smoothing_steps,mesh_parameters,graph,Pnt,maxlev,C);
     corrF=P*corrC;
     corrF(RemoveF)=0;    
     x=x+corrF;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Lev>1: pre-smoothing %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [x,WorkingSet]=ArnoldSmootherContact(WorkingSet,Abc,x,b,Constraint,flip_vector(graph{F},1),mesh_parameters{F}.N_contact,mesh_parameters{F}.Patch_Internal_All,smoothing_steps) ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Computation of the norm of the residual, without considering bcs and active set dofs %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(lev==maxlev)
    resF=b -Abc*x;
    resF(RemoveF)=0;
    resF(find(WorkingSet==1))=0;
    norm_res=norm(resF);
end

end
