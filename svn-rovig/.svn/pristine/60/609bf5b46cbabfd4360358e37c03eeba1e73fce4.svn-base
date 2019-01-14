function [x,WorkingSet,norm_res] = ArnoldVcycleContact2(graph,maxlev,L,maps,x,bnt,mesh,Antnobc,Ant,Pnt, parameters,Constraint,WorkingSet)

EmapGlob2Loc=maps.EmapGlob2Loc;
EmapLoc2Glob=maps.EmapLoc2Glob;
NmapGlob2Loc=maps.NmapGlob2Loc;
NmapLoc2Glob=maps.NmapLoc2Glob;
Patch_Boundary_Edge=maps.Patch_Boundary_Edge;
Patch_Boundary_Node=maps.Patch_Boundary_Node;
Patch_Edge=maps.Patch_Edge;
Patch_Node=maps.Patch_Node;

 % fine and coarse
 F=L;
 C=L-1;
 
 % N and NE on the fine mesh
 NEF=mesh{F}.NE;
 NF=mesh{F}.N;

A=[ Ant{F,1,1}, Ant{F,1,2}, Ant{F,1,3}, Ant{F,1,4};
    Ant{F,2,1}, Ant{F,2,2}, Ant{F,2,3}, Ant{F,2,4};
    Ant{F,3,1}, Ant{F,3,2}, Ant{F,3,3}, Ant{F,3,4};
    Ant{F,4,1}, Ant{F,4,2}, Ant{F,4,3}, Ant{F,4,4}];

Anobc=[ Antnobc{F,1,1}, Antnobc{F,1,2}, Antnobc{F,1,3}, Antnobc{F,1,4};
        Antnobc{F,2,1}, Antnobc{F,2,2}, Antnobc{F,2,3}, Antnobc{F,2,4};
        Antnobc{F,3,1}, Antnobc{F,3,2}, Antnobc{F,3,3}, Antnobc{F,3,4};
        Antnobc{F,4,1}, Antnobc{F,4,2}, Antnobc{F,4,3}, Antnobc{F,4,4}];

if(F==1)
   x=A\bnt; 
else
    


% N and NE on the coarse mesh
 NEC=mesh{C}.NE;
 NC=mesh{C}.N;
 
 %  E_label=mesh{F}.E_label;
%  N_label=mesh{F}.N_label;
 % fine level lenghts
 LFine=[ NEF ; 2 * NEF; 2 * NEF + NF ; 2 * NEF + 2 * NF ];
 LCoarse=[NEC; 2 * NEC;2 * NEC + NC ; 2 * NEC + 2 * NC ];  

 
  E_remove=mesh{F}.E_remove;
  E_dirichlet=mesh{F}.E_dirichlet;

  N_remove=mesh{F}.N_remove;
  N_dirichlet=mesh{F}.N_dirichlet;
  
smoothing_steps=parameters.smoothing_steps;

    
E_remove=mesh{F}.E_remove;
N_remove=mesh{F}.N_remove;    
E_contact_tangentF=mesh{F}.E_contact+NEF;  

removeF=[E_remove,E_remove+NEF,N_remove+2*NEF,N_remove+NF+2*NEF ];
removeF=[removeF,E_contact_tangentF,E_contact_tangentF];
removeF=unique(removeF);
is_on_coarser_grid=0;
top2bottom=1;    



 for ii=1:4
     for jj=1:4
        Antlev{ii,jj}=   Ant{F,ii,jj};
     end
 end
 
energy(1)=0.5*x'*Anobc*x-bnt'*x;
[x,WorkingSet,energy1]=ArnoldSmootherContact6(L,WorkingSet,x,bnt,mesh,Anobc,Antlev,Constraint,smoothing_steps,is_on_coarser_grid,graph,top2bottom,maps);
energy(2)=0.5*x'*Anobc*x-bnt'*x;
% figure
% plot(energy1)

P=Pnt{C};
ActiveSet=find(WorkingSet>0);

% ActiveSet=[];












ActiveSetC=ActiveSet(find(ActiveSet>2*NEF));
ActiveSetC= ActiveSetC(find(ActiveSetC<=2*NEF+NC))-2*NEF+2*NEC;


% no active set
%ActiveSet=[];
P(ActiveSet,:)=0;






% now we want to build coarse gap function
% in particular the constraint for the pressure is always negative
% for the constraint on the displacement in normal direction is take care
% with proper Monotone Restriction
ConstraintF=Constraint-x;
ConstraintC= ArnoldCoarseConstraint(ConstraintF,mesh,maps,C,F);





% compute fine residual on the fine level
resF= bnt - A * x;  

resF(removeF)=0;
%resF(ActiveSet)=0;

resC=P'*resF;
 

E_removeC=mesh{C}.E_remove;
E_removeC1=E_removeC;
E_removeC2=E_removeC1+NEC;

N_removeC=mesh{C}.N_remove;
N_removeC1 = N_removeC + 2 * NEC;
N_removeC2 = N_removeC1+ NC;

E_contact_tangentC=mesh{C}.E_contact+NEC;  

removeC=[E_removeC1,E_removeC2,N_removeC1,N_removeC2];

removeC=[removeC,E_contact_tangentC,E_contact_tangentC];
removeC=unique(removeC);

%removeC=unique([removeC,ActiveSetC']);


%%%%% qui fai Vcycle
corrC=sparse(2*NEC+2*NC,1);
LineSearch=0;


AC=P'*Anobc*P;
AC(removeC,:)=0;
AC(:,removeC)=0;
for ii=removeC
    AC(ii,ii)=1;
end
% AC=0.5*(AC+AC');
resC(removeC,1)=0;

WorkingSet_Loc=[];
B=speye(length(resC));
[corrC,lambda,WorkingSet_Loc] = ArnoldActiveset2(AC,B,resC,ConstraintC,WorkingSet_Loc);


% control if we have energy reduction...

0.5*corrC'*AC*corrC-corrC'*resC;
corrF=P*corrC;

corrF(removeF)=0;


y=x+corrF;

% y=min(y,Constraint);

% c=y-x;
% 
% cAc=c'*Anobc*c;
% 
% % maybe sym(A) in the computation of the residual
% res=bnt-0.5*(A+A')*x;
% res(removeF)=0;
% alpha=min(1,res'*c/cAc);
% x=x+alpha*c;

%good_coarse_corr=(0.5*y'*Anobc*y-bnt'*y<0.5*x'*Anobc*x-bnt'*x)


%0.5*corrF'*Anobc*corrF-resF'*corrF
 x=y;
energy(3)=0.5*x'*Anobc*x-bnt'*x;

is_on_coarser_grid=0;
[x,WorkingSet,energy2]=ArnoldSmootherContact6(L,WorkingSet,x,bnt,mesh,Anobc,Antlev,Constraint,smoothing_steps,is_on_coarser_grid,graph,top2bottom,maps);

energy(4)=0.5*x'*Anobc*x-bnt'*x;

resF=bnt -Anobc*x;
resF(removeF)=0;
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


