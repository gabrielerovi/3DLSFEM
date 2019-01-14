function [x,Nconstraint,Econstraint]=SmootherContact(graph,top2bottom,EmapGlob2Loc,EmapLoc2Glob,NmapGlob2Loc,NmapLoc2Glob,...
    Patch_Boundary_Edge, Patch_Boundary_Node,Patch_Edge, Patch_Node,x,b,mesh,...
    A11_lev,A12_lev,A13_lev,A14_lev,...
    A21_lev,A22_lev,A23_lev,A24_lev,...
    A31_lev,A32_lev,A33_lev,A34_lev,...
    A41_lev,A42_lev,A43_lev,A44_lev,...
    smoothing_steps,is_on_coarser_grid, parameters,Nconstraint,Econstraint)


N=mesh.N;
NE=mesh.NE;
N_remove=mesh.N_remove;
E_remove=mesh.E_remove;
N_removeL=length(N_remove);
E_removeL=length(E_remove);
E_label=mesh.E_label;
N_label=mesh.N_label;
E_bc=mesh.E_bc;
N_bc=mesh.N_bc;
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



% in the 3th column of dirichlet, we have the information on contact 
type_of_dof=1;
[dirichletN,n_and_or_tN,bool_bcN]= boundary_value_bool(type_of_dof);
type_of_dof=2;
[dirichletE,n_and_or_tE,bool_bcE]= boundary_value_bool(type_of_dof);

% we loop on all the edges to find the ones in Gammac
% then we know also the nodes that are in Gammac
contE=0;
contN=0;
for ee=1:NE
    % on boundary
    if(E_bc(ee)>0) 
        % in contact
        if(dirichletE(E_bc(ee),3)==1)
            contE=contE+1;
            Econtact(contE)=ee;
            edge=mesh.edge(ee,:);
            for nn=1:2
            if(dirichletN(N_bc(edge(nn)),3)==1)
            contN=contN+1;
            Ncontact(contN)=edge(nn);
            end
            end
        end
    end
end
    
% we use Arnold smoother on all the vertices that do not belong to GammaC
% so we remove from the vertices the nodes on GammaC
vertices=setdiff(vertices,Ncontact);
    
    
    
%  SMOOTHING-STEPS
for jj=1:smoothing_steps

    % on the vertices we do usual Arnold's smoother
for nn=vertices
% border vertex - edge, Local - Global
    N_Glob=Patch_Node{nn};
    N_BoundaryGlob=Patch_Boundary_Node{nn};
    N_InternalGlob=setdiff(N_Glob,N_BoundaryGlob);
    
%    N_InternalGlob=setdiff(N_InternalGlob,remove_contact_nodes);
    
    E_Glob=Patch_Edge{nn};
    E_BoundaryGlob=Patch_Boundary_Edge{nn};
    E_InternalGlob=setdiff(E_Glob,E_BoundaryGlob);

   % E_InternalGlob=setdiff(E_InternalGlob,remove_contact_edges);
    
    
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
    

   % now we loop on all the edges on GammaC
   
   for ee=Econtact
       % we consider the edge vetices
       edge=mesh.edge(ee,:);
       nn=[];
       % actually only the ones that can be in contact (belong to GammaC)
       for nn1=1:2
       if(dirichletN(N_bc(edge(nn1)),3)==1)
           nn=[nn;edge(nn1)];
       end
       end
    % then we compute Gauss-Seidel for:
    % 1) the displacement of the first node and then check un - g <=0
    % 2) the displacement of the second node and then check un - g <=0
    % 3) the stress of the edge; in frictionless case, then we put (sigma n)_t=0
    %    if both u_i n - g =0, for both node i=1,2, then we put (sigma n ) n = 0 

    
    % Displacement for the k-th node
    for kk=1:2
        
    Tot_InternalGlob=[2 * NE + nn(kk); 2 * NE + N + nn(kk)];
    Tot_ExternalGlob=1:(NE * 2 + N * 2 );
    Tot_ExternalGlob=setdiff(Tot_ExternalGlob,Tot_InternalGlob);
    
    % create local matrix
     A_Loc=[A33_lev(nn(kk),nn(kk)) A34_lev(nn(kk),nn(kk));
            A43_lev(nn(kk),nn(kk)) A44_lev(nn(kk),nn(kk))];   
    % create the rhs 
    b_disp1Loc=  - A((2 * NE + nn(kk)),     Tot_ExternalGlob) * x(Tot_ExternalGlob);
    b_disp2Loc=  - A((2 * NE + N + nn(kk)), Tot_ExternalGlob) * x(Tot_ExternalGlob);
    
    b_disp1Loc= b_disp1Loc+ b_disp1( nn(kk) );
    b_disp2Loc= b_disp2Loc+ b_disp2( nn(kk) ); 
    
    b_Loc=[b_disp1Loc;b_disp2Loc];
    
    % solve the system
    disp=A_Loc\b_Loc;
    
    dof = nn ( kk );
    % update global vector
    [disp,Nconstraint(nn(kk))]=ContactConstraint(mesh,disp,dof,1,parameters);    
    x(Tot_InternalGlob)=disp;
    end
    
    % Stress for the ee edge
    Tot_InternalGlob=[ee; ee + NE];
    Tot_ExternalGlob=1:(NE * 2 + N * 2 );
    Tot_ExternalGlob=setdiff(Tot_ExternalGlob,Tot_InternalGlob);
    % create local matrix
     A_Loc=[A11_lev(ee,ee) A12_lev(ee,ee);
            A21_lev(ee,ee) A22_lev(ee,ee)]; 
        
    % create rhs
    b_sigma1Loc= - A(  ee,              Tot_ExternalGlob) * x(Tot_ExternalGlob);
    b_sigma2Loc= - A(( ee + NE),        Tot_ExternalGlob) * x(Tot_ExternalGlob);

    b_sigma1Loc=b_sigma1Loc+b_sigma1(ee);
    b_sigma2Loc=b_sigma2Loc+b_sigma2(ee);  
    b_Loc=[b_sigma1Loc;b_sigma2Loc];
    
    % solve the system
    sigma=A_Loc\b_Loc;    
    
    % here we enforce:
    % 1) (sigma n)_t = 0  (true for frictionless)
    % 2) (sigma n) n <=0  (always true)

    dof=ee;
    sigma=x_Loc([1,2]);
    % enforce (sigma n) n <=0 and, in frictionless case, (sigma n)_t = 0 
    [sigma,E_sigman_dot_n_less_0]=ContactConstraint(mesh,sigma,dof,type_of_dof,parameters);
        
    % if at least one node does not touch the obstacle, then (n' sigma n)=0    
    if(Nconstraint(edge(1))==0 || Nconstraint(edge(2))==0) 
        Econstraint(ee)=1;
        phidotn=phi_dot_n(mesh,dof);
        normalstress=sigma*phidotn;
        [normalstress_normaldirection,H]=HouseHolderTransformation(normalstress,mesh.normal_edge{dof});
        % since u cdot n = g, n' sigma n= 0
        normalstress_normaldirection(1)=0;    
        % go back to normal stresses
        normalstress=H'*normalstress_normaldirection;   
        % go back to coefficients
        sigma=normalstress./phidotn;
    end

    % update whole solution
    x(Tot_InternalGlob)=sigma;
end

  


end


end

