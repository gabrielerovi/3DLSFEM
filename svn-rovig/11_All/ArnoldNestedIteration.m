function [solF,WorkingSet] = ArnoldNestedIteration(A,Complementarity,h,bnt1,Constraint,mesh_parameters,Pnt)

C=1;
L=length(mesh_parameters);
RemoveNT=mesh_parameters{C}.RemoveNT;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Build the system on the lowest level: A and Complementarity are treated separately   %%%%%%
%%%% Complementarity is a boundary integral, so for consistency on each level:            %%%%%%
%%%% Complementarity_{j}=h_j/h_{j+1} P_{j+1}^j Complementarity_{j+1} P_{j}^{j+1}          %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for lev=L-1:-1:1
A=Pnt{lev}'*A*Pnt{lev};
Complementarity=Pnt{lev}'* Complementarity  *Pnt{lev};
end
Ant1=A+h(1)/h(L)*Complementarity;

Ant1(RemoveNT,:)=0;
for ii=RemoveNT
Ant1(ii,ii)=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Build the coarse constraint starting from the finest one    %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ConstraintC=Constraint;
if(L>1)
for C=L-1:-1:1
F=C+1;
    mp.NEC=mesh_parameters{C}.NE;
    mp.NEF=mesh_parameters{F}.NE;
    mp.NC=mesh_parameters{C}.N;
    mp.NF=mesh_parameters{F}.N;
    mp.EC_contact=mesh_parameters{C}.E_contact;
    mp.Patch_Edge_Monotone=mesh_parameters{C}.Patch_Edge_Monotone;
    mp.Patch_Node_Monotone=mesh_parameters{C}.Patch_Node_Monotone;   
ConstraintC= ArnoldCoarseConstraint(ConstraintC,mp);
end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Solve the coarsest system, with the proper projected rhs    %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WorkingSet_Loc=[];
B=speye(length(bnt1));
[solC,lambda,WorkingSet] = ArnoldActiveset(Ant1,B,bnt1,ConstraintC,WorkingSet_Loc);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Interpolate the coarse workingset and the solution          %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for lev=1:L-1
C=lev;
F=lev+1;
solC=Pnt{lev}*solC;
    mp.NEC=mesh_parameters{C}.NE;
    mp.NEF=mesh_parameters{F}.NE;
    mp.NC=mesh_parameters{C}.N;
    mp.NF=mesh_parameters{F}.N;
    mp.EC_contact=mesh_parameters{C}.E_contact;
    mp.Patch_Edge_Monotone=mesh_parameters{C}.Patch_Edge_Monotone;
    mp.Patch_Node_Monotone=mesh_parameters{C}.Patch_Node_Monotone;
WorkingSet= WSC2WSF(WorkingSet,mp);
end

solF=solC;
solF=min(solF,Constraint);


end