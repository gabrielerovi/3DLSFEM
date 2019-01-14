function [sol_k,is_the_right_solution,Constraint]=activeset4(Ant,b,Constraint,mesh,actual_solution,mapping_loc)
 
toll=10^(-15);
E_contact=mesh.E_contact;
N_contact=mesh.N_contact;
% E_remove=mesh.E_remove;
% N_remove=mesh.N_remove;
% E_label=mesh.E_label;
% N_label=mesh.N_label;
NE=mesh.NE;
N=mesh.N;

CheckConstraintE1=Constraint.CheckConstraintE1;
RhsE1=Constraint.RhsE1;
CheckConstraintN1=Constraint.CheckConstraintN1;
RhsN1=Constraint.RhsN1;
CheckConstraintE2=Constraint.CheckConstraintE2;
RhsE2=Constraint.RhsE2;
CheckConstraintN2=Constraint.CheckConstraintN2;
RhsN2=Constraint.RhsN2;


WorkingSetE=Constraint.WorkingSetE;
WorkingSetN=Constraint.WorkingSetN;

NotWorkingSetE=setdiff(E_contact,WorkingSetE);
NotWorkingSetN=setdiff(N_contact,WorkingSetN);

    

 
%  for ii=1:2
%  
%      for jj=1:4
%           Ant{ii,jj}(E_remove,:)=0;
%      end
% 
%      for jj=1:4
%      Ant{ii+2,jj}(N_remove,:)=0;
%      end
%  
%  
%      for jj=E_remove
%           Ant{ii,ii}(jj,jj)=1;
%      end
%  
%      for jj=N_remove
%      Ant{ii+2,ii+2}(jj,jj)=1;
%      end
%  
%  end
%  
% 
%  
%  EcontBC=0;
%  type_of_dof=2;
%  for ii=E_remove
%  U_xy=[0;0];
%  EcontBC=EcontBC+1;
%  tmp=add_boundary_bc_elasticity2D(U_xy,type_of_dof,E_label(EcontBC), 0);
%  coeff = RT_dirichlet_coeff(E_remove(EcontBC), mesh);
%  jj1=ii;
%  jj2=ii+NE;
%  b(jj1)=tmp(1)/coeff;
%  b(jj2)=tmp(2)/coeff;
%  end
%  
%  NcontBC=0;
%  type_of_dof=1;
%  for ii=N_remove
%  U_xy=[0;0];
%  NcontBC=NcontBC+1;
%  tmp=add_boundary_bc_elasticity2D(U_xy,type_of_dof,N_label(NcontBC), 0);
%  jj1=ii+2*NE;
%  jj2=ii+2*NE+N;
%  b(jj1)=tmp(1);
%  b(jj2)=tmp(2);
%  end
 
 
 
 
 
 
 
 
 
 

 
 
 
 
 
 
 
 
 

% INPUT: 1) the 16 blocks of the matrix Ant without boundary conditions
%           Ant is A where, on GammaC, we have normal and tangent components
%           instead of the canonical ones
%        2) WorkingSetE: on these edges we enforce sigma_n =0
%        3) WorkingSetN: on these nodes we enforce u_n=0
%
% we solve for the system built in this way and then we check:
% if (CheckConstraintE<RhsE && CheckConstraintN<RhsN)
% then we have found the solution
% otherwise consider the dofs related to:
% 1) WorkingSetE:
% if the midpoint penetrates into the obstacle: then, remove it from WorkingSetE
% 2) E_contact - WorkingSetE:
% if the dof is positive, add it to WorkingSetE
% 3) WorkingSetN:
% if the stress in the vertex is bigger than zero: then remove it from WorkingSetN
% 4) N_contact - WorkingSetN:
% if the dof is bigger than the gap, add it to WorkingSetN

 is_the_right_solution=false;
WorkingSetE_loc_numbering=[];
for ii=1:length(WorkingSetE)
WorkingSetE_loc_numbering(ii)=find(WorkingSetE(ii)==E_contact);
end
NotWorkingSetE_loc_numbering=setdiff(1:length(E_contact),WorkingSetE_loc_numbering);

WorkingSetN_loc_numbering=[];
for ii=1:length(WorkingSetN)
WorkingSetN_loc_numbering(ii)=find(WorkingSetN(ii)==N_contact);
end
NotWorkingSetN_loc_numbering=setdiff(1:length(N_contact),WorkingSetN_loc_numbering);


% ENFORCE FRICTIONLESS 
% the 2-2 block is the tangent one for stresses
if(~isempty(E_contact))
for ee=E_contact
    
    for ii=1:4
    Ant{2,ii}(ee,:)=0; 
    end
    
Ant{2,2}(ee,ee)=1;
b(NE+ee,1)=0;
end
end
% ENFORCE ZERO NORMAL STRESS ON WorkingSetE
% the 1-1 block is the normal one for stresses
contWE=0;
if(~isempty(WorkingSetE))
for ee=WorkingSetE
    for ii=1:4
    Ant{1,ii}(ee,:)=0;
    end
contWE=contWE+1;   
Ant{1,1}(ee,ee)=1;
b(ee,1)=RhsE2(WorkingSetE_loc_numbering(contWE));
end
end
% ENFORCE DISPLACEMENT==GAP ON WorkingSetN
% the 3-3 block is the normal one for displacements
contWN=0;
if(~isempty(WorkingSetN))
for nn=WorkingSetN
    for ii=1:4
    Ant{3,ii}(nn,:)=0;
    end
contWN=contWN+1;
Ant{3,3}(nn,nn)=1;
b(2*NE+nn,1)=RhsN2(WorkingSetN_loc_numbering(contWN));
end
end

A_k=[Ant{1,1} Ant{1,2} Ant{1,3}  Ant{1,4};
     Ant{2,1} Ant{2,2} Ant{2,3}  Ant{2,4};
     Ant{3,1} Ant{3,2} Ant{3,3}  Ant{3,4};
     Ant{4,1} Ant{4,2} Ant{4,3}  Ant{4,4}    ];
 
sol_k_tmp=A_k\b;
 
actual_solution(mapping_loc)=sol_k_tmp;
sol_k=actual_solution;

% REMOVE AND ADD EDGE CONSTRAINTS
% in WorkingSetE/N and in E/N_contact we have an enumeration wrt all the
% edges/nodes, while in CheckConstraintE/N we have only the dofs in contact
% so we have to take care of such a local numberings
RemoveFromConstraintE=find(CheckConstraintE1(WorkingSetE_loc_numbering,:) * sol_k - RhsE1(WorkingSetE_loc_numbering)>toll);
RemoveFromConstraintE=unique(E_contact(WorkingSetE_loc_numbering(RemoveFromConstraintE)));

% RemoveFromConstraintE2=find(CheckConstraintE2(WorkingSetE_loc_numbering,:) * sol_k - RhsE2(WorkingSetE_loc_numbering)>toll);
% RemoveFromConstraintE=[RemoveFromConstraintE1;RemoveFromConstraintE2];

% AddToConstraintE1=find(CheckConstraintE1(NotWorkingSetE_loc_numbering,:) * sol_k - RhsE1(NotWorkingSetE_loc_numbering)>toll);
AddToConstraintE=find(CheckConstraintE2(NotWorkingSetE_loc_numbering,:) * sol_k - RhsE2(NotWorkingSetE_loc_numbering)>toll);
AddToConstraintE=unique(E_contact(NotWorkingSetE_loc_numbering(AddToConstraintE)));
AddToConstraintE=setdiff(AddToConstraintE,WorkingSetE);
% AddToConstraintE=[AddToConstraintE1;AddToConstraintE2];

% RemoveFromConstraintE=intersect(RemoveFromConstraintE,WorkingSetE);

% AddToConstraintE=setdiff(AddToConstraintE,WorkingSetE);
% 

WorkingSetE=setdiff(WorkingSetE,RemoveFromConstraintE);
Constraint.WorkingSetE=[WorkingSetE,AddToConstraintE];


% REMOVE AND ADD NODE CONSTRAINTS


RemoveFromConstraintN=find(CheckConstraintN1(WorkingSetN_loc_numbering,:) * sol_k - RhsN1(WorkingSetN_loc_numbering)>toll);
RemoveFromConstraintN=unique(N_contact(WorkingSetN_loc_numbering(RemoveFromConstraintN)));

% RemoveFromConstraintN2=find(CheckConstraintN2(WorkingSetN_loc_numbering,:) * sol_k - RhsN2(WorkingSetN_loc_numbering)>toll);
% RemoveFromConstraintN=[RemoveFromConstraintN1;RemoveFromConstraintN2];

% AddToConstraintN1=find(CheckConstraintN1(NotWorkingSetN_loc_numbering,:) * sol_k - RhsN1(NotWorkingSetN_loc_numbering)>toll);
AddToConstraintN=find(CheckConstraintN2(NotWorkingSetN_loc_numbering,:) * sol_k - RhsN2(NotWorkingSetN_loc_numbering)>toll);
AddToConstraintN=unique(N_contact(NotWorkingSetN_loc_numbering(AddToConstraintN)));
AddToConstraintN=setdiff(AddToConstraintN,WorkingSetN);

% AddToConstraintN=[AddToConstraintN1;AddToConstraintN2];


% RemoveFromConstraintN=intersect(RemoveFromConstraintN,WorkingSetN);
 
% AddToConstraintN=setdiff(AddToConstraintN,WorkingSetN);

WorkingSetN=setdiff(WorkingSetN,RemoveFromConstraintN);
Constraint.WorkingSetN=[WorkingSetN,AddToConstraintN];


% if we are not changing the working set, it means all the constraints are satisfied
% the enforced bc were correct and we are done
if(isempty(RemoveFromConstraintE) && isempty(RemoveFromConstraintN) && isempty(AddToConstraintE) && isempty(AddToConstraintN) )
 is_the_right_solution=true;
end
 
sol_k=sol_k_tmp;
end
 
