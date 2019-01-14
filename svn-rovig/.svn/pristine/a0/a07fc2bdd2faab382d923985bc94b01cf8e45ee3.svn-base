
function x=ArnoldSmoother7(A,x,b,Patch_Internal_All,smoothing_steps,graph,top2bottom) 


if(top2bottom)
    vertices=graph;
else
    vertices=fliplr(graph);
end

%  SMOOTHING-STEPS
for jj=1:smoothing_steps

    xold=x;
    
 % on interior nodes we use standard Arnold patch smoother
for nn=vertices

    Patch=Patch_Internal_All{nn};
  
    A_Loc= A (Patch,Patch);
    b_Loc= b (Patch)- A(Patch,:) * x;
    c_Loc=A_Loc\b_Loc;

    x(Patch)=x(Patch) + c_Loc;
 
end

end


end





















