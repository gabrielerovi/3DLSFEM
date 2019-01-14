


function x=ArnoldSmoother4(graph,top2bottom,EmapGlob2Loc,EmapLoc2Glob,NmapGlob2Loc,NmapLoc2Glob,...
    Patch_Boundary_Edge, Patch_Boundary_Node,Patch_Edge, Patch_Node,x,b,mesh,...
    A11_lev,A12_lev,A21_lev,A22_lev,...
    smoothing_steps,is_on_coarser_grid)

N=mesh.N;
NE=mesh.NE;
N_remove=mesh.N_remove;
E_remove=mesh.E_remove;
N_removeL=length(N_remove);
E_removeL=length(E_remove);
E_label=mesh.E_label;
N_label=mesh.N_label;

LGlob=[NE; 2*NE];


A=[ A11_lev A12_lev;
    A21_lev A22_lev];

b_sigma1=b(1:LGlob(1));
b_sigma2=b(1+LGlob(1):LGlob(2));




if(is_on_coarser_grid==true)
for ii=E_remove
    b_sigma1(ii)=0;
    b_sigma2(ii)=0;
end 
end



if(top2bottom)
    vertices=graph;
else
    vertices=fliplr(graph);
end



%  SMOOTHING-STEPS
for jj=1:smoothing_steps

for nn=vertices

    % border vertex - edge, Local - Global
    N_Glob=Patch_Node{nn};
    N_BoundaryGlob=Patch_Boundary_Node{nn};
    N_InternalGlob=setdiff(N_Glob,N_BoundaryGlob);
    
    E_Glob=Patch_Edge{nn};
    E_BoundaryGlob=Patch_Boundary_Edge{nn};
    E_InternalGlob=setdiff(E_Glob,E_BoundaryGlob);

    
    % consider ALL THE EDGES
    %E_InternalGlob=E_Glob;
    
    
    NELoc=length(E_Glob);
    NLoc= length(N_Glob);
    
    
    % all vertex -edge, Local - Global
    N_Loc=1:NLoc;
    N_BoundaryLoc=cell2mat(values(NmapGlob2Loc{nn},num2cell(N_BoundaryGlob,1)));
    N_InternalLoc=cell2mat(values(NmapGlob2Loc{nn},num2cell(N_InternalGlob,1)));  

    E_Loc=1:NELoc;   
    E_BoundaryLoc=cell2mat(values(EmapGlob2Loc{nn},num2cell(E_BoundaryGlob,1)));
    E_InternalLoc=cell2mat(values(EmapGlob2Loc{nn},num2cell(E_InternalGlob,1)));


    LLoc=[NELoc; NELoc+NLoc];
    
    a11=A11_lev(E_InternalGlob,E_InternalGlob); 
    a12=A12_lev(E_InternalGlob,E_InternalGlob);
    
    a21=A21_lev(E_InternalGlob,E_InternalGlob); 
    a22=A22_lev(E_InternalGlob,E_InternalGlob); 


    E_ExternalGlob=1:NE;
    E_ExternalGlob=setdiff(E_ExternalGlob,E_InternalGlob);
    E_ExternalGlob1=E_ExternalGlob;
    E_ExternalGlob2=E_ExternalGlob+NE;
    E_InternalGlob1=E_InternalGlob;
    E_InternalGlob2=E_InternalGlob1+NE;
    
    Tot_ExternalGlob=[E_ExternalGlob1,E_ExternalGlob2];
    Tot_InternalGlob=[E_InternalGlob1,E_InternalGlob2];
    % assign the global external forces to the rhs
    b_sigma1Loc= - A(E_InternalGlob1,     Tot_ExternalGlob) * x(Tot_ExternalGlob);
    b_sigma1Loc=b_sigma1Loc+b_sigma1(E_InternalGlob);

    b_sigma2Loc= - A(E_InternalGlob2,     Tot_ExternalGlob) * x(Tot_ExternalGlob);
    b_sigma2Loc=b_sigma2Loc+b_sigma2(E_InternalGlob);
    
    
    A_Loc=[a11 a12; 
           a21 a22];
       
    b_Loc=[b_sigma1Loc;b_sigma2Loc];
    
    x_Loc=A_Loc\b_Loc;
    
    xold=x;
    x(Tot_InternalGlob)=x_Loc;
%     xold-x
%     norm(xold-x)
    
    
end



end


end

