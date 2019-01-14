
function [x,WorkingSet]=ArnoldSmootherContactAugmented(WorkingSet,A,x,b,Constraint,vertices,vertices_GammaC,Patch_Internal_All,smoothing_steps) 


vertices_interior=setdiff(vertices,vertices_GammaC);
clearvars vertices

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
    
    b_Loc= b (Patch)- A(Patch,:) * x;
    c_Loc=A_Loc\b_Loc;
    x(Patch)=x(Patch) + c_Loc;
    

end


for nn=vertices_GammaC

    Patch=Patch_Internal_All{nn};
  
    A_Loc= A (Patch,Patch);
    b_Loc= b (Patch)- A(Patch,:) * x;
 
    constraint_Loc=Constraint(Patch)-x(Patch);
         
     B_Loc=speye(length(A_Loc));
     
     WorkingSet_Loc=find(WorkingSet(Patch)>0);

     
    [c_Loc,lambda,WorkingSet_Loc] =  ActiveSet(A_Loc,B_Loc,b_Loc,constraint_Loc,WorkingSet_Loc,zeros(length(b_Loc),1))  ;   
    WorkingSet (Patch(setdiff(1:length(b_Loc),WorkingSet_Loc) ))=0;
    WorkingSet (Patch(WorkingSet_Loc) )=1;
     
      x(Patch)=x(Patch) + c_Loc;
%  if(min(constraint_Loc-c_Loc)<0)
%      jj,nn
%      fermamiplease=1
%  end
end


end


end





















