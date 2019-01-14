function [solF,WorkingSet] = ArnoldNestedIteration3D(A,Complementarity,h,bnt1,Constraint,mesh_parameters,Pnt)

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
Ant1=A+(h(1)/h(L))^2*Complementarity;

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
    mp.NFC=mesh_parameters{C}.NF;
    mp.NFF=mesh_parameters{F}.NF;
    mp.NC=mesh_parameters{C}.N;
    mp.NF=mesh_parameters{F}.N;
    mp.FC_contact=mesh_parameters{C}.F_contact;
    mp.Patch_Face_Monotone=mesh_parameters{C}.Patch_Face_Monotone;
    mp.Patch_Node_Monotone=mesh_parameters{C}.Patch_Node_Monotone;   
ConstraintC= ArnoldCoarseConstraint3D(ConstraintC,mp);
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
    mp.NFC=mesh_parameters{C}.NF;
    mp.NFF=mesh_parameters{F}.NF;
    mp.NC=mesh_parameters{C}.N;
    mp.NF=mesh_parameters{F}.N;
    mp.FC_contact=mesh_parameters{C}.F_contact;
    mp.Patch_Face_Monotone=mesh_parameters{C}.Patch_Face_Monotone;
    mp.Patch_Node_Monotone=mesh_parameters{C}.Patch_Node_Monotone;
WorkingSet= WSC2WSF3D(WorkingSet,mp);
end

solF=solC;
solF=min(solF,Constraint);


end