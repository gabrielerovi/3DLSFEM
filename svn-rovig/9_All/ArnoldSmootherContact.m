
function [x,WorkingSet]=ArnoldSmootherContact(WorkingSet,A,x,b,Constraint,vertices,vertices_GammaC,Patch_Internal_All,smoothing_steps) 


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
%     c_Loc=A_Loc\b_Loc;

[L,U] = lu(full(A (Patch,Patch)));
yy=L\b_Loc;
c_Loc=U\yy;

for hh=1:1
res=b_Loc-A_Loc*c_Loc;
yy=L\res;
cc=U\yy;
c_Loc=c_Loc+cc;
end

x(Patch)=x(Patch) + c_Loc;


end


for nn=vertices_GammaC

    Patch=Patch_Internal_All{nn};
  
    A_Loc= A (Patch,Patch);
    b_Loc= b (Patch)- A(Patch,:) * x;
    c_Loc=A_Loc\b_Loc;
 
    constraint_Loc=Constraint(Patch)-x(Patch);
         
     B_Loc=speye(length(A_Loc));
     
     WorkingSet_Loc=find(WorkingSet(Patch)>0);

     [c_Loc,lambda,WorkingSet_Loc] = ArnoldActiveset(A_Loc,B_Loc,b_Loc,constraint_Loc,WorkingSet_Loc);
          
    WorkingSet (Patch(setdiff(1:length(b_Loc),WorkingSet_Loc) ))=0;
    WorkingSet (Patch(WorkingSet_Loc) )=1;
     
    x(Patch)=x(Patch) + c_Loc;
 
end

end


end





















