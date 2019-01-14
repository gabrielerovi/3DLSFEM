


function x=ArnoldSmoother2(graph,top2bottom,EmapGlob2Loc,EmapLoc2Glob,NmapGlob2Loc,NmapLoc2Glob,...
    Patch_Boundary_Edge, Patch_Boundary_Node,Patch_Edge, Patch_Node,x,b,mesh,...
    A11_lev,smoothing_steps,eta,is_on_coarser_grid)

N=mesh.N;
NE=mesh.NE;
N_remove=mesh.N_remove;
E_remove=mesh.E_remove;
N_removeL=length(N_remove);
E_removeL=length(E_remove);
E_label=mesh.E_label;
N_label=mesh.N_label;

LGlob=[NE;2*NE;2*NE+N;2*NE+2*N ];


b_sigma1=b(1:LGlob(1));


if(is_on_coarser_grid==true)
for ii=E_remove
    b_sigma1(ii)=0;
end
   
end



if(top2bottom)
    vertices=graph;
else
    vertices=fliplr(graph);
end



%  SMOOTHING-STEPS
for jj=1:smoothing_steps

x_old=x;    
sigma1=x(1:LGlob(1));

 for ii=E_remove
    sigma1(ii)=b_sigma1(ii);
 end

for nn=vertices

    % border vertex - edge, Local - Global
    N_Glob=Patch_Node{nn};
    N_BoundaryGlob=Patch_Boundary_Node{nn};
    N_InternalGlob=setdiff(N_Glob,N_BoundaryGlob);
    
    E_Glob=Patch_Edge{nn};
    E_BoundaryGlob=Patch_Boundary_Edge{nn};
    E_InternalGlob=setdiff(E_Glob,E_BoundaryGlob);

    NELoc=length(E_Glob);
    NLoc= length(N_Glob);
    
    
    % all vertex -edge, Local - Global
    N_Loc=1:NLoc;
    N_BoundaryLoc=cell2mat(values(NmapGlob2Loc{nn},num2cell(N_BoundaryGlob,1)));
    N_InternalLoc=cell2mat(values(NmapGlob2Loc{nn},num2cell(N_InternalGlob,1)));  

    E_Loc=1:NELoc;   
    E_BoundaryLoc=cell2mat(values(EmapGlob2Loc{nn},num2cell(E_BoundaryGlob,1)));
    E_InternalLoc=cell2mat(values(EmapGlob2Loc{nn},num2cell(E_InternalGlob,1)));


    LLoc=[NELoc; 2*NELoc; 2*NELoc+NLoc;2*NELoc+2*NLoc;];
    
    a11=A11_lev(E_Glob,E_Glob); 
 
    % assign the global external forces to the rhs
    b_sigma1Loc=b_sigma1(E_Glob);
    
    b_sigma1Loc(E_BoundaryLoc)=sigma1(E_BoundaryGlob);
    
    A_Loc=[a11];
            
    diagBC=[E_BoundaryLoc];
    A_Loc(diagBC,:)=0;
    for iii=diagBC
        A_Loc(iii,iii)=1;
    end
    
    b_Loc=[b_sigma1Loc];
    
    for bbb=diagBC
        
    end

    x_Loc=A_Loc\b_Loc;
    sigma1Loc=x_Loc(1:LLoc(1));

    sigma1(E_InternalGlob)= sigma1Loc(E_InternalLoc); 

end

x=[sigma1];

norm(x_old-x);

end


end

