
function [x,WorkingSet,energy]=ArnoldSmootherContact7(WorkingSet,Aenergy,A,x,b,Constraint,mesh,Patch_Internal_All,smoothing_steps,graph,top2bottom) 


if(top2bottom)
    vertices=graph;
else
    vertices=fliplr(graph);
end
% interior nodes that do not lie on GammaC
vertices_GammaC=mesh.N_contact;
vertices_interior=setdiff(vertices,vertices_GammaC);


contEnergy=0;
%  SMOOTHING-STEPS
for jj=1:smoothing_steps

    xold=x;
    
 % on interior nodes we use standard Arnold patch smoother
for nn=vertices_interior

    Patch=Patch_Internal_All{nn};
  
    A_Loc= A (Patch,Patch);
    b_Loc= b (Patch)- A(Patch,:) * x;
    c_Loc=A_Loc\b_Loc;

    x(Patch)=x(Patch) + c_Loc;
    
    contEnergy=contEnergy+1;
    energy(contEnergy)=0.5*x'*Aenergy*x-b'*x;

end


for nn=vertices_GammaC

    Patch=Patch_Internal_All{nn};
  
    A_Loc= A (Patch,Patch);
    b_Loc= b (Patch)- A(Patch,:) * x;
    c_Loc=A_Loc\b_Loc;
 
    constraint_Loc=Constraint(Patch)-x(Patch);
         
     B_Loc=speye(length(A_Loc));
     
     WorkingSet_Loc=find(WorkingSet(Patch)>0);

     [c_Loc,lambda,WorkingSet_Loc] = ArnoldActiveset2(A_Loc,B_Loc,b_Loc,constraint_Loc,WorkingSet_Loc);
     
%      [c_Loc,WorkingSet_Loc]=activesetResidual2(A_Loc,b_Loc,constraint_Loc,WorkingSet_Loc);
     
    WorkingSet (Patch(setdiff(1:length(b_Loc),WorkingSet_Loc) ))=0;
    WorkingSet (Patch(WorkingSet_Loc) )=1;
     
    x(Patch)=x(Patch) + c_Loc;
    contEnergy=contEnergy+1;
    energy(contEnergy)=0.5*x'*Aenergy*x-b'*x;    
 
end

end


end





















