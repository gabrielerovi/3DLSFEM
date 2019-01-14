


function [x,WorkingSetNormal_N,WorkingSetNormal_E]=ArnoldSmootherContactNormalTangent(graph,top2bottom,EmapGlob2Loc,EmapLoc2Glob,NmapGlob2Loc,NmapLoc2Glob,...
    Patch_Boundary_Edge, Patch_Boundary_Node,Patch_Edge, Patch_Node,x,b,mesh,...
    A11_lev,A12_lev,A13_lev,A14_lev,...
    A21_lev,A22_lev,A23_lev,A24_lev,...
    A31_lev,A32_lev,A33_lev,A34_lev,...
    A41_lev,A42_lev,A43_lev,A44_lev,...
    smoothing_steps,is_on_coarser_grid)

N=mesh.N;
NE=mesh.NE;

WorkingSetNormal_N=sparse(N,1);
WorkingSetNormal_E=sparse(NE,1);


N_remove=mesh.N_remove;
E_remove=mesh.E_remove;
N_removeL=length(N_remove);
E_removeL=length(E_remove);
E_label=mesh.E_label;
N_label=mesh.N_label;

LGlob=[NE;2*NE; 2*NE+N; 2*NE+ 2*N ];


A=[A11_lev A12_lev A13_lev A14_lev;
   A21_lev A22_lev A23_lev A24_lev;
   A31_lev A32_lev A33_lev A34_lev;
   A41_lev A42_lev A43_lev A44_lev;
    ];

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
% interior nodes that do not lie on GammaC
vertices_GammaC=mesh.N_contact;
vertices_interior=setdiff(vertices,vertices_GammaC);


%  SMOOTHING-STEPS
for jj=1:smoothing_steps

    
 % on interior nodes we use standard Arnold patch smoother
for nn=vertices_interior

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
    
    a11=A11_lev(E_InternalGlob,E_InternalGlob); 
    a12=A12_lev(E_InternalGlob,E_InternalGlob); 
    a13=A13_lev(E_InternalGlob,N_InternalGlob); 
    a14=A14_lev(E_InternalGlob,N_InternalGlob); 
    
    a21=A21_lev(E_InternalGlob,E_InternalGlob); 
    a22=A22_lev(E_InternalGlob,E_InternalGlob); 
    a23=A23_lev(E_InternalGlob,N_InternalGlob); 
    a24=A24_lev(E_InternalGlob,N_InternalGlob); 
    
    a31=A31_lev(N_InternalGlob,E_InternalGlob); 
    a32=A32_lev(N_InternalGlob,E_InternalGlob); 
    a33=A33_lev(N_InternalGlob,N_InternalGlob); 
    a34=A34_lev(N_InternalGlob,N_InternalGlob);
    
    a41=A41_lev(N_InternalGlob,E_InternalGlob); 
    a42=A42_lev(N_InternalGlob,E_InternalGlob); 
    a43=A43_lev(N_InternalGlob,N_InternalGlob); 
    a44=A44_lev(N_InternalGlob,N_InternalGlob); 
    


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
    
    
    
    A_Loc=[a11 a12 a13 a14;
           a21 a22 a23 a24;
           a31 a32 a33 a34;
           a41 a42 a43 a44;];
    b_Loc=[b_sigma1Loc;b_sigma2Loc;b_disp1Loc;b_disp2Loc];
    x_Loc=A_Loc\b_Loc;
    x(Tot_InternalGlob)=x_Loc;

end










% on interior nodes we use standard Arnold patch smoother
for nn=vertices_GammaC

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
    
    a11=A11_lev(E_InternalGlob,E_InternalGlob); 
    a12=A12_lev(E_InternalGlob,E_InternalGlob); 
    a13=A13_lev(E_InternalGlob,N_InternalGlob); 
    a14=A14_lev(E_InternalGlob,N_InternalGlob); 
    
    a21=A21_lev(E_InternalGlob,E_InternalGlob); 
    a22=A22_lev(E_InternalGlob,E_InternalGlob); 
    a23=A23_lev(E_InternalGlob,N_InternalGlob); 
    a24=A24_lev(E_InternalGlob,N_InternalGlob); 
    
    a31=A31_lev(N_InternalGlob,E_InternalGlob); 
    a32=A32_lev(N_InternalGlob,E_InternalGlob); 
    a33=A33_lev(N_InternalGlob,N_InternalGlob); 
    a34=A34_lev(N_InternalGlob,N_InternalGlob);
    
    a41=A41_lev(N_InternalGlob,E_InternalGlob); 
    a42=A42_lev(N_InternalGlob,E_InternalGlob); 
    a43=A43_lev(N_InternalGlob,N_InternalGlob); 
    a44=A44_lev(N_InternalGlob,N_InternalGlob); 
    


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
    
    
    
    A_loc=[a11 a12 a13 a14;
           a21 a22 a23 a24;
           a31 a32 a33 a34;
           a41 a42 a43 a44;];
    b_loc=[b_sigma1Loc;b_sigma2Loc;b_disp1Loc;b_disp2Loc];
    
    NE_loc=length(b_sigma1Loc);
    N_loc=length(b_disp1Loc);
    Tot_loc=2*NE_loc+2*N_loc;
    
    ST_loc=sparse(0,Tot_loc);
    CT_loc=sparse(0,0);
    BN_loc=sparse(0,Tot_loc);
    CN_loc=sparse(0,0);
    
    % by construction the node nn is on GammaC
    N_contact=nn;
    E_contact=mesh.E_bc_contact(E_InternalGlob1);
    
    BN_loc(1,2*NE_loc+1)=-1;
    CN_loc(1,1)=-gap_function(mesh.node(nn,1), mesh.node(nn,2));
        
    contE=0;
    cont=0;
    for ee=1:NE_loc
        cont=cont+1;
        if(E_contact(cont)>0)
            contE=contE+1;
            E_really_in_contact(contE)=E_contact(cont);
            ST_loc(contE, NE_loc+cont)=1;
            CT_loc(contE,1)=0;
            % the first one is for the normal displacement
            phidotn=phi_dot_n(mesh,E_InternalGlob(ee));
            BN_loc(1+contE,cont)=-1*sign(phidotn);
            CN_loc(1+contE,1)=0;
        end        
    end
    
    [x_Loc,lambda,WorkingSet_tmp_loc] = activeset(A_loc,ST_loc',BN_loc',b_loc,ST_loc,CT_loc,BN_loc,CN_loc,x(Tot_InternalGlob));
    

    WorkingSet_loc=sparse(length(E_really_in_contact)+1,1);
    WorkingSet_loc(WorkingSet_tmp_loc)=1;
    
    WorkingSetNormal_N(nn)=WorkingSet_loc(1);
    
    WorkingSetNormal_E(E_really_in_contact)=WorkingSet_loc(2:end);
        
    
    
    x(Tot_InternalGlob)=x_Loc;

end







end


end

