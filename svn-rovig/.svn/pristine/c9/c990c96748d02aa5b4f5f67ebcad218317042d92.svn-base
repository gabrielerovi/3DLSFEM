
function [x,WorkingSet,energy]=ArnoldSmootherContact6(L,WorkingSet,x,b,mesh,Aenergy,Ant,Constraint,smoothing_steps,is_on_coarser_grid,graph,top2bottom,maps)

grid=mesh;
mesh=mesh{L};
householder=1;
[M_Normal_Tangent,M_Normal_TangentT] = MatrixOnGammaCwithNormalTangentComponents(mesh,householder);

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
E_contact=mesh.E_contact;
N_contact=mesh.N_contact;


big_value=10^10;

% c=big_value*ones(2*NE+2*N,1);
% 
% 
% c(E_contact)=Constraint.RhsE2;
% c(N_contact+2*NE)=Constraint.RhsN2;

LGlob=[NE;2*NE; 2*NE+N; 2*NE+ 2*N ];


    A=    [Ant{1,1} Ant{1,2} Ant{1,3} Ant{1,4};
           Ant{2,1} Ant{2,2} Ant{2,3} Ant{2,4};
           Ant{3,1} Ant{3,2} Ant{3,3} Ant{3,4};
           Ant{4,1} Ant{4,2} Ant{4,3} Ant{4,4};];
    
    Asym=0.5 * (A+A');
    
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

contEnergy=0;
%  SMOOTHING-STEPS
for jj=1:smoothing_steps

    xold=x;
    
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
%     b_sigma1Loc= - Asym(E_InternalGlob1,     :) * x(:);
%     b_sigma2Loc= - Asym(E_InternalGlob2,     :) * x(:);
%     b_disp1Loc=  - Asym(N_InternalGlob1,     :) * x(:);
%     b_disp2Loc=  - Asym(N_InternalGlob2,     :) * x(:);
    b_sigma1Loc= - A(E_InternalGlob1,     :) * x(:);
    b_sigma2Loc= - A(E_InternalGlob2,     :) * x(:);
    b_disp1Loc=  - A(N_InternalGlob1,     :) * x(:);
    b_disp2Loc=  - A(N_InternalGlob2,     :) * x(:);
    
    
    
    b_sigma1Loc=b_sigma1Loc+b_sigma1(E_InternalGlob);
    b_sigma2Loc=b_sigma2Loc+b_sigma2(E_InternalGlob);
    b_disp1Loc=b_disp1Loc+b_disp1(N_InternalGlob);
    b_disp2Loc=b_disp2Loc+b_disp2(N_InternalGlob);    
    
    
    
    A_tot_Loc=[A_Loc{1,1} A_Loc{1,2} A_Loc{1,3} A_Loc{1,4};
           A_Loc{2,1} A_Loc{2,2} A_Loc{2,3} A_Loc{2,4};
           A_Loc{3,1} A_Loc{3,2} A_Loc{3,3} A_Loc{3,4};
           A_Loc{4,1} A_Loc{4,2} A_Loc{4,3} A_Loc{4,4};];
       
    b_Loc=[b_sigma1Loc;b_sigma2Loc;b_disp1Loc;b_disp2Loc];
    c_Loc=A_tot_Loc\b_Loc;
    x(Tot_InternalGlob)=x(Tot_InternalGlob) + c_Loc;
   
    contEnergy=contEnergy+1;
    energy(contEnergy)=0.5*x'*Aenergy*x-b'*x;

end




for nn=vertices_GammaC
    % border vertex - edge, Local - Global
    N_Glob=Patch_Node{nn};
    N_BoundaryGlob=Patch_Boundary_Node{nn};
    N_InternalGlob=setdiff(N_Glob,N_BoundaryGlob);
    
    E_Glob=Patch_Edge{nn};
    E_BoundaryGlob=Patch_Boundary_Edge{nn};
    E_InternalGlob=setdiff(E_Glob,E_BoundaryGlob);
    
    N_contact_Glob=intersect(N_contact,N_InternalGlob); 
    E_contact_Glob=intersect(E_contact,E_InternalGlob); 
    
%     % since the friction dof is zero, we
%     E_InternalGlob2=setdiff(E_InternalGlob,E_contact_Glob)
    
    
%     cont=0;
%     E_contact_Loc=[];
%     for ii=E_contact_Glob
%         cont=cont+1;
%         E_contact_Loc(cont)=find(E_InternalGlob==ii);
%     end
%     cont=0;
%     N_contact_Loc=[];
%     for ii=N_contact_Glob
%         cont=cont+1;
%         N_contact_Loc(cont)=find(N_InternalGlob==ii);
%     end 
    
    NELoc=length(E_InternalGlob);
    NLoc= length(N_InternalGlob);
    
    % contact only in normal component (tangent is frictionless for edges)
%     Contact_Loc=[E_contact_Loc, N_contact_Loc+2*NELoc];
    
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
%     b_sigma1Loc= - Asym(E_InternalGlob1,     :) * x(:);
%     b_sigma2Loc= - Asym(E_InternalGlob2,     :) * x(:);
%     b_disp1Loc=  - Asym(N_InternalGlob1,     :) * x(:);
%     b_disp2Loc=  - Asym(N_InternalGlob2,     :) * x(:);
    
    b_sigma1Loc= - A(E_InternalGlob1,     :) * x(:);
    b_sigma2Loc= - A(E_InternalGlob2,     :) * x(:);
    b_disp1Loc=  - A(N_InternalGlob1,     :) * x(:);
    b_disp2Loc=  - A(N_InternalGlob2,     :) * x(:);
    
    
    b_sigma1Loc=b_sigma1Loc+b_sigma1(E_InternalGlob);
    b_sigma2Loc=b_sigma2Loc+b_sigma2(E_InternalGlob);
    b_disp1Loc=b_disp1Loc+b_disp1(N_InternalGlob);
    b_disp2Loc=b_disp2Loc+b_disp2(N_InternalGlob);    
    
    
    
    A_tot_Loc=[A_Loc{1,1} A_Loc{1,2} A_Loc{1,3} A_Loc{1,4};
           A_Loc{2,1} A_Loc{2,2} A_Loc{2,3} A_Loc{2,4};
           A_Loc{3,1} A_Loc{3,2} A_Loc{3,3} A_Loc{3,4};
           A_Loc{4,1} A_Loc{4,2} A_Loc{4,3} A_Loc{4,4};];
       
    b_Loc=[b_sigma1Loc;b_sigma2Loc;b_disp1Loc;b_disp2Loc];
    
    c=Constraint(Tot_InternalGlob)-x(Tot_InternalGlob);
    
    %[jj,nn]
    %  c_Loc = ArnoldMinimumEnergy (A_tot_Loc,b_Loc,c, Contact_Loc);
     
     B=speye(length(A_tot_Loc));
     
     WorkingSet_Loc=find(WorkingSet(Tot_InternalGlob)>0);
     [c_Loc,lambda,WorkingSet_Loc] = ArnoldActiveset2(A_tot_Loc,B,b_Loc,c,WorkingSet_Loc);
     
    WorkingSet (Tot_InternalGlob(setdiff(1:2*NELoc+2*NLoc,WorkingSet_Loc) ))=0;
    WorkingSet (Tot_InternalGlob(WorkingSet_Loc) )=1;
     
     x(Tot_InternalGlob)=x(Tot_InternalGlob) + c_Loc;
   
    contEnergy=contEnergy+1;
    energy(contEnergy)=0.5*x'*Aenergy*x-b'*x;
    
%     sol=M_Normal_Tangent*x;
%     print_displacement_solution(grid,sol(1+2*grid{L}.NE:2*grid{L}.NE+grid{L}.N)',sol(1+2*grid{L}.NE+grid{L}.N:end)');


end







end
% figure
% plot(energy,'ro');


end





































% 
% 
% 
% function [x,Constraint]=ArnoldSmootherContact6(AFinenobcnt,AArnoldLocal,L,activeset,x,b,mesh,Ant,Constraint,smoothing_steps,is_on_coarser_grid,graph,top2bottom,maps)
% 
% meshprova=mesh;
% mesh=mesh{L};
% 
% householder=1;
% [M_Normal_Tangent,M_Normal_TangentT] = MatrixOnGammaCwithNormalTangentComponents(mesh,householder);
% 
% 
% graph=graph{L};
% N_dirichlet=mesh.N_dirichlet;
% E_dirichlet=mesh.E_dirichlet;
% EmapGlob2Loc=maps.EmapGlob2Loc{L};
% EmapLoc2Glob=maps.EmapLoc2Glob{L};
% NmapGlob2Loc=maps.NmapGlob2Loc{L};
% NmapLoc2Glob=maps.NmapLoc2Glob{L};
% Patch_Boundary_Edge=maps.Patch_Boundary_Edge{L};
% Patch_Boundary_Node=maps.Patch_Boundary_Node{L};
% Patch_Edge=maps.Patch_Edge{L};
% Patch_Node=maps.Patch_Node{L};
% 
% N=mesh.N;
% NE=mesh.NE;
% N_remove=mesh.N_remove;
% E_remove=mesh.E_remove;
% E_contact=mesh.E_contact;
% N_contact=mesh.N_contact;
% 
% 
% big_value=10^10;
% 
% c=big_value*ones(2*NE+2*N,1);
% 
% 
% c(E_contact)=Constraint.RhsE2;
% c(N_contact+2*NE)=Constraint.RhsN2;
% 
% LGlob=[NE;2*NE; 2*NE+N; 2*NE+ 2*N ];
% 
% 
%     A=    [Ant{1,1} Ant{1,2} Ant{1,3} Ant{1,4};
%            Ant{2,1} Ant{2,2} Ant{2,3} Ant{2,4};
%            Ant{3,1} Ant{3,2} Ant{3,3} Ant{3,4};
%            Ant{4,1} Ant{4,2} Ant{4,3} Ant{4,4};];
%     
%     Asym=0.5 * (A+A');
%     
% b_sigma1=b(1:LGlob(1));
% b_sigma2=b(1+LGlob(1):LGlob(2));
% b_disp1=b(1+LGlob(2):LGlob(3));
% b_disp2=b(1+LGlob(3):LGlob(4));
% 
% if(is_on_coarser_grid==true)
%     
% for ii=E_remove
%     b_sigma1(ii)=0;
%     b_sigma2(ii)=0;
% end
% 
% for ii=N_remove
%     b_disp1(ii)=0;
%     b_disp2(ii)=0;
% end
% end
% 
% 
% if(top2bottom)
%     vertices=graph;
% else
%     vertices=fliplr(graph);
% end
% % interior nodes that do not lie on GammaC
% vertices_GammaC=mesh.N_contact;
% vertices_interior=setdiff(vertices,vertices_GammaC);
% 
% 
% %  SMOOTHING-STEPS
% for jj=1:smoothing_steps
% 
%     xold=x;
%       
%  
% for nn=1:N
% nn
%     % border vertex - edge, Local - Global
%     N_Glob=Patch_Node{nn};
%     N_BoundaryGlob=Patch_Boundary_Node{nn};
%     
%     % we consider all the nodes that can be in contact so that, when we
%     % apply the constraints, they will be considered
%     
%     N_InternalGlob=N_Glob;
% 
%    
%     E_Glob=Patch_Edge{nn};
%     E_BoundaryGlob=Patch_Boundary_Edge{nn};
%     E_InternalGlob=E_Glob;
%     
%     
%     
%     NELoc=length(E_Glob);
%     NLoc= length(N_Glob);
%     
%     
%     % all vertex -edge, Local - Global
%     N_Loc=1:NLoc;
%     N_BoundaryLoc=cell2mat(values(NmapGlob2Loc{nn},num2cell(N_BoundaryGlob,1)));
%     N_InternalLoc=cell2mat(values(NmapGlob2Loc{nn},num2cell(N_InternalGlob,1)));  
% 
%     E_Loc=1:NELoc;   
%     E_BoundaryLoc=cell2mat(values(EmapGlob2Loc{nn},num2cell(E_BoundaryGlob,1)));
%     E_InternalLoc=cell2mat(values(EmapGlob2Loc{nn},num2cell(E_InternalGlob,1)));
% 
% 
%     LLoc=[NELoc; 2 * NELoc; 2 * NELoc + NLoc; 2 * NELoc + 2 * NLoc];
%     
% 
%     
% 
%     E_ExternalGlob=1:NE;
%     E_ExternalGlob=setdiff(E_ExternalGlob,E_InternalGlob);
%     E_ExternalGlob1=E_ExternalGlob;
%     E_ExternalGlob2=E_ExternalGlob+NE;
%     E_InternalGlob1=E_InternalGlob;
%     E_InternalGlob2=E_InternalGlob1+NE;    
%     
%     N_ExternalGlob=1:N;
%     N_ExternalGlob=setdiff(N_ExternalGlob,N_InternalGlob);
%     N_ExternalGlob1=N_ExternalGlob+2*NE;
%     N_ExternalGlob2=N_ExternalGlob1+N;
%     
%     N_InternalGlob1=N_InternalGlob+2*NE;
%     N_InternalGlob2=N_InternalGlob1+N; 
%     
%     N_BoundaryGlob1=N_BoundaryGlob+2*NE;
%     N_BoundaryGlob2=N_BoundaryGlob1+N;
%         
%     
%     Tot_ExternalGlob=[E_ExternalGlob1,E_ExternalGlob2,N_ExternalGlob1,N_ExternalGlob2];
%     Tot_InternalGlob=[E_InternalGlob1,E_InternalGlob2,N_InternalGlob1,N_InternalGlob2];
%     
%     b_sigma1Loc= b_sigma1(E_InternalGlob) ;
%     b_sigma2Loc= b_sigma2(E_InternalGlob) ;
%     b_disp1Loc= b_disp1(N_InternalGlob) ;
%     b_disp2Loc= b_disp2(N_InternalGlob) ;  
% 
%    
%     
%     cont=0;
%     for ii=N_BoundaryLoc
%         cont=cont+1;
%         b_disp1Loc(ii)=x(N_BoundaryGlob1(cont));
%         b_disp2Loc(ii)=x(N_BoundaryGlob2(cont));        
%     end
%     
%     b_Loc=[b_sigma1Loc;b_sigma2Loc;b_disp1Loc;b_disp2Loc];
% 
%     % enforce bc on the local matrix...
%     N_bc=find(N_dirichlet(N_Glob)==1);
%     E_bc=find(E_dirichlet(E_Glob)==1);
%     remove=[E_bc,E_bc+NELoc, N_bc+2*NELoc,N_bc+2*NELoc+NLoc];
%     
%     a_Loc=AArnoldLocal{nn};
%     
%     a_Loc(remove,:)=0;
%     for rr=remove
%         a_Loc(rr,rr)=1;
%     end
%   
%     c_Loc=c(Tot_InternalGlob);
%     
%     tmpE=intersect(activeset(find(activeset<=NE)),E_InternalGlob);
%     tmpN=intersect(activeset(find(activeset>NE)),N_InternalGlob);
%     
%     activesetE_Loc= cell2mat(values(EmapGlob2Loc{nn},num2cell(tmpE)));
%     activesetN_Loc= cell2mat(values(NmapGlob2Loc{nn},num2cell(tmpN)));
%     
%     activeset_Loc=[activesetE_Loc, activesetN_Loc+2*NELoc];
%   
%     nn,jj
%     activeset
%     
%     activeset=setdiff(activeset, activeset_Loc);
%     [x_Loc,activeset_Loc]=activesetResidual(a_Loc,b_Loc,c_Loc,activeset_Loc);
% 
%     
%     activeset=unique([activeset,Tot_InternalGlob(activeset_Loc)]);
%     
%     
%     y=x;
% %     ED=find(E_dirichlet==1);
% %     ND=find(N_dirichlet==1);
% %     rem=[ED,ED+ NE,  ND + 2 * NE,  ND + 2 * NE + N];
% %     
% %     y(rem)=0;
%     energy_1=0.5* y'*AFinenobcnt *y - b'*y;
% 
%     x(Tot_InternalGlob')=x_Loc;
% 
% %     sol=M_Normal_Tangent*x;
% %     print_displacement_solution(meshprova,sol(1+2*meshprova{L}.NE:2*meshprova{L}.NE+meshprova{L}.N)',sol(1+2*meshprova{L}.NE+meshprova{L}.N:end)');
% 
% 
%     y=x;
%     energy_2=0.5* y'*AFinenobcnt *y - b'*y;
%     
%     diff=energy_1-energy_2
%     if(diff<0)
%         fermami=1
%     end
%     energy(nn)=energy_2;
% end
% 
% 
% norm(x-xold)
% 
% end
% 
% 
% 
% end

























































% Tot_ExternalGlob=1:2*NE+2*N;
% for nn=N_contact
%     
%     tmp1=[nn,nn+N]+2*NE;
%     tmp=setdiff(Tot_ExternalGlob,tmp1);
%     
%     
%     b_disp(1,1)= b_disp1(nn)  - A(nn+2*NE,tmp) *x(tmp) ;
%     b_disp(2,1)= b_disp2(nn)  - A(nn+2*NE+N,tmp) *x(tmp) ;
%     
%     A_disp=[Ant{3,3}(nn,nn), Ant{3,3}(nn,nn); 
%             Ant{4,3}(nn,nn), Ant{4,4}(nn,nn);];
%    
%     x_Loc=A_disp\b_disp;
%     
%     c_Loc=x_Loc-x(tmp1);
% %     
% %     cAc=c_Loc'*A_disp*c_Loc;
% %     
% %     res=  
% %     alpha= c_Loc'*()
%     if(x_Loc(1)>c(nn+2*NE))
%         A_disp(1,1)=1;
%         A_disp(1,2)=0;
%         b_disp(1)=c(nn+2*NE);
%         x_Loc=A_disp\b_disp;
%     end
%     
%     x(tmp1)=x_Loc;
%     
% end
% 
% 
% 
% for ee=E_contact
%     
%     tmp1=[ee,ee+NE];
%     tmp=setdiff(Tot_ExternalGlob,tmp1);
%     
%     b_sigma(1,1)= b_sigma1(ee)  - A(ee,tmp) *x(tmp) ;
%     b_sigma(2,1)= b_sigma2(ee)  - A(ee+NE,tmp) *x(tmp) ;
%     
%     A_sigma=[Ant{1,1}(ee,ee), Ant{1,2}(ee,ee); 
%              Ant{2,1}(ee,ee), Ant{2,2}(ee,ee);];
%    
%     x_Loc=A_sigma\b_sigma;
%     
% 
%     if(x_Loc(1)>c(ee))
%         A_sigma(1,1)=1;
%         b_sigma(1,2)=0;
%         b_sigma(1)=c(ee);
%         x_Loc=A_sigma\b_sigma;
%     end
%     
%     x(tmp1)=x_Loc;
%     
% end























%%% WE SOLVE FOR INTERNAL AND THEN ON CONTACT
% NE=mesh.NE;
% N=mesh.N;
% Ec=mesh.E_contact;
% Nc=mesh.N_contact;
% dofs=[Ec,Ec+NE,Nc+2*NE,Nc+2*NE+N];
% dofs_ext=1:2*NE+2*N;
% dofs_ext=setdiff(dofs_ext,dofs);
% 
% b_sigma1Loc= b_sigma1(Ec) - A(Ec,     dofs_ext) * x(dofs_ext) ;
% b_sigma2Loc= b_sigma2(Ec) - A(Ec,     dofs_ext) * x(dofs_ext) ;
% b_disp1Loc= b_disp1(Nc)   - A(Nc,     dofs_ext) * x(dofs_ext) ;
% b_disp2Loc= b_disp2(Nc)   - A(Nc,     dofs_ext) * x(dofs_ext) ;
% b_Loc=[b_sigma1Loc;b_sigma2Loc;b_disp1Loc;b_disp2Loc];
% 
% a_Loc=[ Ant{1,1}(Ec,Ec), Ant{1,2}(Ec,Ec), Ant{1,3}(Ec,Nc), Ant{1,4}(Ec,Nc);
%         Ant{2,1}(Ec,Ec), Ant{2,2}(Ec,Ec), Ant{2,3}(Ec,Nc), Ant{2,4}(Ec,Nc);
%         Ant{3,1}(Nc,Ec), Ant{3,2}(Nc,Ec), Ant{3,3}(Nc,Nc), Ant{3,4}(Nc,Nc);
%         Ant{4,1}(Nc,Ec), Ant{4,2}(Nc,Ec), Ant{4,3}(Nc,Nc), Ant{4,4}(Nc,Nc);];
%     
% 
% c_Loc=c(dofs);
% 
% tmpE=find(activeset<=NE);
% tmpN=find(activeset>NE);
%     
% 
% 
% activeset_Loc=[tmpE,tmpN];
% [x_Loc,activeset_Loc]=activesetResidual(a_Loc,b_Loc,c_Loc,activeset_Loc);
% 
% 
% x(dofs)=x_Loc;





% % on interior nodes we use standard Arnold patch smoother
% for nn=vertices_GammaC
% 
%     % border vertex - edge, Local - Global
%     N_Glob=Patch_Node{nn};
%     N_BoundaryGlob=Patch_Boundary_Node{nn};
%     
%     % we consider all the nodes that can be in contact so that, when we
%     % apply the constraints, they will be considered
%     
%     % N_InternalGlob=nn;
%     N_InternalGlob=nn;%N_Glob;
%    
%     E_Glob=Patch_Edge{nn};
%     E_BoundaryGlob=Patch_Boundary_Edge{nn};
%     E_InternalGlob=setdiff(E_Glob,E_BoundaryGlob);
%     E_InternalGlob=E_Glob;
%     
%     
%     
%     NELoc=length(E_Glob);
%     NLoc= length(N_Glob);
%     
%     
%     % all vertex -edge, Local - Global
%     N_Loc=1:NLoc;
%     N_BoundaryLoc=cell2mat(values(NmapGlob2Loc{nn},num2cell(N_BoundaryGlob,1)));
%     N_InternalLoc=cell2mat(values(NmapGlob2Loc{nn},num2cell(N_InternalGlob,1)));  
% 
%     E_Loc=1:NELoc;   
%     E_BoundaryLoc=cell2mat(values(EmapGlob2Loc{nn},num2cell(E_BoundaryGlob,1)));
%     E_InternalLoc=cell2mat(values(EmapGlob2Loc{nn},num2cell(E_InternalGlob,1)));
% 
% 
%     LLoc=[NELoc; 2 * NELoc; 2 * NELoc + NLoc; 2 * NELoc + 2 * NLoc];
%     
%     A_Loc{1,1}=Ant{1,1}(E_InternalGlob,E_InternalGlob); 
%     A_Loc{1,2}=Ant{1,2}(E_InternalGlob,E_InternalGlob); 
%     A_Loc{1,3}=Ant{1,3}(E_InternalGlob,N_InternalGlob); 
%     A_Loc{1,4}=Ant{1,4}(E_InternalGlob,N_InternalGlob); 
%     
%     A_Loc{2,1}=Ant{2,1}(E_InternalGlob,E_InternalGlob); 
%     A_Loc{2,2}=Ant{2,2}(E_InternalGlob,E_InternalGlob); 
%     A_Loc{2,3}=Ant{2,3}(E_InternalGlob,N_InternalGlob); 
%     A_Loc{2,4}=Ant{2,4}(E_InternalGlob,N_InternalGlob); 
%     
%     A_Loc{3,1}=Ant{3,1}(N_InternalGlob,E_InternalGlob); 
%     A_Loc{3,2}=Ant{3,2}(N_InternalGlob,E_InternalGlob); 
%     A_Loc{3,3}=Ant{3,3}(N_InternalGlob,N_InternalGlob); 
%     A_Loc{3,4}=Ant{3,4}(N_InternalGlob,N_InternalGlob);
%     
%     A_Loc{4,1}=Ant{4,1}(N_InternalGlob,E_InternalGlob); 
%     A_Loc{4,2}=Ant{4,2}(N_InternalGlob,E_InternalGlob); 
%     A_Loc{4,3}=Ant{4,3}(N_InternalGlob,N_InternalGlob); 
%     A_Loc{4,4}=Ant{4,4}(N_InternalGlob,N_InternalGlob); 
%     
% 
% 
%     E_ExternalGlob=1:NE;
%     E_ExternalGlob=setdiff(E_ExternalGlob,E_InternalGlob);
%     E_ExternalGlob1=E_ExternalGlob;
%     E_ExternalGlob2=E_ExternalGlob+NE;
%     E_InternalGlob1=E_InternalGlob;
%     E_InternalGlob2=E_InternalGlob1+NE;    
%     
%     N_ExternalGlob=1:N;
%     N_ExternalGlob=setdiff(N_ExternalGlob,N_InternalGlob);
%     N_ExternalGlob1=N_ExternalGlob+2*NE;
%     N_ExternalGlob2=N_ExternalGlob1+N;
%     N_InternalGlob1=N_InternalGlob+2*NE;
%     N_InternalGlob2=N_InternalGlob1+N; 
%     
%     Tot_ExternalGlob=[E_ExternalGlob1,E_ExternalGlob2,N_ExternalGlob1,N_ExternalGlob2];
%     Tot_InternalGlob=[E_InternalGlob1,E_InternalGlob2,N_InternalGlob1,N_InternalGlob2];
% 
%     b_sigma1Loc= b_sigma1(E_InternalGlob)- A(E_InternalGlob1,     :) *x ;
%     b_sigma2Loc= b_sigma2(E_InternalGlob)- A(E_InternalGlob2,     :) *x ;
%     b_disp1Loc= b_disp1(N_InternalGlob)  - A(N_InternalGlob1,     :) *x ;
%     b_disp2Loc= b_disp2(N_InternalGlob)  - A(N_InternalGlob2,     :) *x ;  
%     
% %     b_sigma1Loc= b_sigma1(E_InternalGlob)- Asym(E_InternalGlob1,     Tot_ExternalGlob) *x(Tot_ExternalGlob) ;
% %     b_sigma2Loc= b_sigma2(E_InternalGlob)- Asym(E_InternalGlob2,     Tot_ExternalGlob) *x(Tot_ExternalGlob) ;
% %     b_disp1Loc= b_disp1(N_InternalGlob)  - Asym(N_InternalGlob1,     Tot_ExternalGlob) *x(Tot_ExternalGlob) ;
% %     b_disp2Loc= b_disp2(N_InternalGlob)  - Asym(N_InternalGlob2,     Tot_ExternalGlob) *x(Tot_ExternalGlob) ;
%     
%     a_Loc=[];
%     for ii=1:4
%     a_Loc=[a_Loc;A_Loc{ii,1}, A_Loc{ii,2}, A_Loc{ii,3}, A_Loc{ii,4};];
%     end
%     
%     
%     b_Loc=[b_sigma1Loc;b_sigma2Loc;b_disp1Loc;b_disp2Loc];
%   
%      c_Loc=c(Tot_InternalGlob)-x(Tot_InternalGlob);
% %    c_Loc=c(Tot_InternalGlob);
%     
%     tmpE=intersect(activeset(find(activeset<=NE)),E_InternalGlob);
%     tmpN=intersect(activeset(find(activeset>NE)),N_InternalGlob);
%     
%     activesetE_Loc= cell2mat(values(EmapGlob2Loc{nn},num2cell(tmpE)));
%     activesetN_Loc= cell2mat(values(NmapGlob2Loc{nn},num2cell(tmpN)));
%     
%     activeset_Loc=[activesetE_Loc, activesetN_Loc+2*NELoc];
%   
%     nn,jj
%     activeset
%     
%     activeset=setdiff(activeset, Tot_InternalGlob);
%     [x_Loc,activeset_Loc]=activesetResidual(a_Loc,b_Loc,c_Loc,activeset_Loc);
% 
%     activeset=unique([activeset,Tot_InternalGlob(activeset_Loc)]);
%     
%     x(Tot_InternalGlob')=x(Tot_InternalGlob')+x_Loc;
% 
% end
% 
% 
% 
% 
% 
% norm(x-xold)
% 
% end




