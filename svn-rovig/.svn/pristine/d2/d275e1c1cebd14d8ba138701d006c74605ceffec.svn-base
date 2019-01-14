
function [x,WorkingSet]=ArnoldSmootherContactEdge(WorkingSet,A,x,b,Constraint,edges_interior,edges_GammaC,Edge_Patch_Internal_All,smoothing_steps) 


lengthEInterior=length(edges_interior);
lengthEGammaC=length(edges_GammaC);
%  SMOOTHING-STEPS
for jj=1:smoothing_steps

    xold=x;
 % on interior nodes we use standard Arnold patch smoother

for ee1=1:lengthEInterior

    ee=edges_interior(ee1);
    Patch=Edge_Patch_Internal_All{ee};
  
    A_Loc= A (Patch,Patch);
    b_Loc= b (Patch)- A(Patch,:) * x;
    c_Loc=A_Loc\b_Loc;

    x(Patch)=x(Patch) + c_Loc;
    
    b_Loc= b (Patch)- A(Patch,:) * x;
    c_Loc=A_Loc\b_Loc;
    x(Patch)=x(Patch) + c_Loc;
    

end



for ee2=1:lengthEGammaC

    ee=edges_GammaC(ee2);
    Patch=Edge_Patch_Internal_All{ee};
  
    A_Loc= A (Patch,Patch);
    b_Loc= b (Patch)- A(Patch,:) * x;
 
    constraint_Loc=Constraint(Patch)-x(Patch);
         
     B_Loc=speye(length(A_Loc));
     
     WorkingSet_Loc=find(WorkingSet(Patch)>0);

     
    [c_Loc,lambda,WorkingSet_Loc] =  ActiveSet(A_Loc,B_Loc,b_Loc,constraint_Loc,WorkingSet_Loc,zeros(length(b_Loc),1))  ;   
    WorkingSet (Patch(setdiff(1:length(b_Loc),WorkingSet_Loc) ))=0;
    WorkingSet (Patch(WorkingSet_Loc) )=1;
     
      x(Patch)=x(Patch) + c_Loc;

end


end


end





















