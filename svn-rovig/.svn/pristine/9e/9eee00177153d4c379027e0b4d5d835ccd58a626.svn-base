function [solF,WorkingSet] = ArnoldNestedIteration(Ant,bnt_lev,Constraint,mesh,maps,Pnt)

C=1;
L=length(mesh);
RemoveNT=mesh{C}.RemoveNT;

% MATRIX
AC=[ Ant{C,1,1}, Ant{C,1,2}, Ant{C,1,3}, Ant{C,1,4};
    Ant{C,2,1}, Ant{C,2,2}, Ant{C,2,3}, Ant{C,2,4};
    Ant{C,3,1}, Ant{C,3,2}, Ant{C,3,3}, Ant{C,3,4};
    Ant{C,4,1}, Ant{C,4,2}, Ant{C,4,3}, Ant{C,4,4}];

AC(RemoveNT,:)=0;
for ii=RemoveNT
AC(ii,ii)=1;
end

% RHS
bC=bnt_lev{C};

% CONSTRAINT
ConstraintC=Constraint;
if(L>1)
for C=L-1:-1:1
F=C+1;
ConstraintC= ArnoldCoarseConstraint(ConstraintC,mesh,maps,C,F);
end    
end


WorkingSet_Loc=[];
B=speye(length(bC));
[solC,lambda,WorkingSet] = ArnoldActiveset2(AC,B,bC,ConstraintC,WorkingSet_Loc);
% [solC2,WorkingSetC2]=activesetResidual2(AC,bC,ConstraintC,WorkingSet_Loc);


% [M_Normal_Tangent,M_Normal_TangentT] = MatrixOnGammaCwithNormalTangentComponents(mesh{C},1);
% sol=M_Normal_Tangent*solC;
% print_displacement_solution2(C,mesh,sol(1+2*mesh{C}.NE:2*mesh{C}.NE+mesh{C}.N)',sol(1+2*mesh{C}.NE+mesh{C}.N:end)');
% sol=M_Normal_Tangent*solC2;
% print_displacement_solution2(C,mesh,sol(1+2*mesh{C}.NE:2*mesh{C}.NE+mesh{C}.N)',sol(1+2*mesh{C}.NE+mesh{C}.N:end)');

for lev=1:L-1
C=lev;
F=lev+1;
solC=Pnt{lev}*solC;
WorkingSet= WSC2WSF(WorkingSet,mesh,maps,C,F);
end

solF=solC;
solF=min(solF,Constraint);
% [M_Normal_Tangent,M_Normal_TangentT] = MatrixOnGammaCwithNormalTangentComponents(mesh{C},1);
% sol=M_Normal_Tangent*solC;
% print_displacement_solution2(C,mesh,sol(1+2*mesh{C}.NE:2*mesh{C}.NE+mesh{C}.N)',sol(1+2*mesh{C}.NE+mesh{C}.N:end)');

% [M_Normal_Tangent,M_Normal_TangentT] = MatrixOnGammaCwithNormalTangentComponents(mesh{L},1);
% sol=M_Normal_Tangent*solF;
% print_displacement_solution2(L,mesh,sol(1+2*mesh{L}.NE:2*mesh{L}.NE+mesh{L}.N)',sol(1+2*mesh{L}.NE+mesh{L}.N:end)');


end