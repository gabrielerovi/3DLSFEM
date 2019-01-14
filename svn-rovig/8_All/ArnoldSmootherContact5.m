


function [x]=ArnoldSmootherContact5(L,x,b,mesh,Ant,smoothing_steps,is_on_coarser_grid,graph,top2bottom,maps)




meshprova=mesh;


[M_Normal_Tangent,M_Normal_TangentT] = MatrixOnGammaCwithNormalTangentComponents(mesh,1);



graph=graph{L};
EmapGlob2Loc=maps.EmapGlob2Loc{L};
EmapLoc2Glob=maps.EmapLoc2Glob{L};
NmapGlob2Loc=maps.NmapGlob2Loc{L};
NmapLoc2Glob=maps.NmapLoc2Glob{L};
Patch_Boundary_Edge=maps.Patch_Boundary_Edge{L};
Patch_Boundary_Node=maps.Patch_Boundary_Node{L};
Patch_Edge=maps.Patch_Edge{L};
Patch_Node=maps.Patch_Node{L};




N=mesh.N;
NE=mesh.NE;
N_remove=mesh.N_remove;
E_remove=mesh.E_remove;
N_removeL=length(N_remove);
E_removeL=length(E_remove);
E_label=mesh.E_label;
N_label=mesh.N_label;

LGlob=[NE;2*NE; 2*NE+N; 2*NE+ 2*N ];


    A=    [Ant{1,1} Ant{1,2} Ant{1,3} Ant{1,4};
           Ant{2,1} Ant{2,2} Ant{2,3} Ant{2,4};
           Ant{3,1} Ant{3,2} Ant{3,3} Ant{3,4};
           Ant{4,1} Ant{4,2} Ant{4,3} Ant{4,4};];
       
b_sigma1=b(1:LGlob(1));
b_sigma2=b(1+LGlob(1):LGlob(2));
b_disp1=b(1+LGlob(2):LGlob(3));
b_disp2=b(1+LGlob(3):LGlob(4));

if(is_on_coarser_grid==true)
    
for ii=E_remove
    b_sigma1(ii)=0;
    b_sigma2(ii)=0;
end

for ii=N_remove
    b_disp1(ii)=0;
    b_disp2(ii)=0;
end
end


if(top2bottom)
    vertices=graph;
else
    vertices=fliplr(graph);
end


%  SMOOTHING-STEPS
for jj=1:smoothing_steps

    
 % on interior nodes we use standard Arnold patch smoother
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


    LLoc=[NELoc; 2 * NELoc; 2 * NELoc + NLoc; 2 * NELoc + 2 * NLoc];
    
    A_Loc{1,1}=Ant{1,1}(E_InternalGlob,E_InternalGlob); 
    A_Loc{1,2}=Ant{1,2}(E_InternalGlob,E_InternalGlob); 
    A_Loc{1,3}=Ant{1,3}(E_InternalGlob,N_InternalGlob); 
    A_Loc{1,4}=Ant{1,4}(E_InternalGlob,N_InternalGlob); 
    
    A_Loc{2,1}=Ant{2,1}(E_InternalGlob,E_InternalGlob); 
    A_Loc{2,2}=Ant{2,2}(E_InternalGlob,E_InternalGlob); 
    A_Loc{2,3}=Ant{2,3}(E_InternalGlob,N_InternalGlob); 
    A_Loc{2,4}=Ant{2,4}(E_InternalGlob,N_InternalGlob); 
    
    A_Loc{3,1}=Ant{3,1}(N_InternalGlob,E_InternalGlob); 
    A_Loc{3,2}=Ant{3,2}(N_InternalGlob,E_InternalGlob); 
    A_Loc{3,3}=Ant{3,3}(N_InternalGlob,N_InternalGlob); 
    A_Loc{3,4}=Ant{3,4}(N_InternalGlob,N_InternalGlob);
    
    A_Loc{4,1}=Ant{4,1}(N_InternalGlob,E_InternalGlob); 
    A_Loc{4,2}=Ant{4,2}(N_InternalGlob,E_InternalGlob); 
    A_Loc{4,3}=Ant{4,3}(N_InternalGlob,N_InternalGlob); 
    A_Loc{4,4}=Ant{4,4}(N_InternalGlob,N_InternalGlob); 
    


    E_ExternalGlob=1:NE;
    E_ExternalGlob=setdiff(E_ExternalGlob,E_InternalGlob);
    E_ExternalGlob1=E_ExternalGlob;
    E_ExternalGlob2=E_ExternalGlob+NE;
    E_InternalGlob1=E_InternalGlob;
    E_InternalGlob2=E_InternalGlob1+NE;    
    
    N_ExternalGlob=1:N;
    N_ExternalGlob=setdiff(N_ExternalGlob,N_InternalGlob);
    N_ExternalGlob1=N_ExternalGlob+2*NE;
    N_ExternalGlob2=N_ExternalGlob1+N;
    N_InternalGlob1=N_InternalGlob+2*NE;
    N_InternalGlob2=N_InternalGlob1+N; 
    
    
    Tot_ExternalGlob=[E_ExternalGlob1,E_ExternalGlob2,N_ExternalGlob1,N_ExternalGlob2];
    Tot_InternalGlob=[E_InternalGlob1,E_InternalGlob2,N_InternalGlob1,N_InternalGlob2];
    % assign the global external forces to the rhs
    b_sigma1Loc= - A(E_InternalGlob1,     Tot_ExternalGlob) * x(Tot_ExternalGlob);
    b_sigma2Loc= - A(E_InternalGlob2,     Tot_ExternalGlob) * x(Tot_ExternalGlob);
    b_disp1Loc=  - A(N_InternalGlob1,     Tot_ExternalGlob) * x(Tot_ExternalGlob);
    b_disp2Loc=  - A(N_InternalGlob2,     Tot_ExternalGlob) * x(Tot_ExternalGlob);
     
    b_sigma1Loc=b_sigma1Loc+b_sigma1(E_InternalGlob);
    b_sigma2Loc=b_sigma2Loc+b_sigma2(E_InternalGlob);
    b_disp1Loc=b_disp1Loc+b_disp1(N_InternalGlob);
    b_disp2Loc=b_disp2Loc+b_disp2(N_InternalGlob);    
    
    
    
    A_tot_Loc=[A_Loc{1,1} A_Loc{1,2} A_Loc{1,3} A_Loc{1,4};
           A_Loc{2,1} A_Loc{2,2} A_Loc{2,3} A_Loc{2,4};
           A_Loc{3,1} A_Loc{3,2} A_Loc{3,3} A_Loc{3,4};
           A_Loc{4,1} A_Loc{4,2} A_Loc{4,3} A_Loc{4,4};];
       
    b_Loc=[b_sigma1Loc;b_sigma2Loc;b_disp1Loc;b_disp2Loc];
    x_Loc=A_tot_Loc\b_Loc;
    x(Tot_InternalGlob)=x_Loc;
    
%     sol=M_Normal_Tangent*x;
%     print_displacement_solution(meshprova,sol(1+2*meshprova{L}.NE:2*meshprova{L}.NE+meshprova{L}.N)',sol(1+2*meshprova{L}.NE+meshprova{L}.N:end)');
end

end


end

