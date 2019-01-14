


function [x,Constraint]=ArnoldSmootherContact3(x,b,mesh,Ant,Constraint,smoothing_steps,is_on_coarser_grid,graph,top2bottom,maps)



L=length(mesh);

meshprova=mesh;

mesh=mesh{L};

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
% interior nodes that do not lie on GammaC
vertices_GammaC=mesh.N_contact;
vertices_interior=setdiff(vertices,vertices_GammaC);

vertices_GammaC=[vertices_GammaC,flip(vertices_GammaC)];
vertices_interior=[vertices_interior,flip(vertices_interior)];

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

%     sol=M_Normal_Tangent*x;
%     print_displacement_solution(meshprova,sol(1+2*meshprova{L}.NE:2*meshprova{L}.NE+meshprova{L}.N)',sol(1+2*meshprova{L}.NE+meshprova{L}.N:end)');









% on interior nodes we use standard Arnold patch smoother
for nn=vertices_GammaC

    % border vertex - edge, Local - Global
    N_Glob=Patch_Node{nn};
    N_BoundaryGlob=Patch_Boundary_Node{nn};
    
    % we consider all the nodes that can be in contact so that, when we
    % apply the constraints, they will be considered
    
    % N_InternalGlob=nn;
    N_InternalGlob=nn;%N_Glob;
   
    E_Glob=Patch_Edge{nn};
    E_BoundaryGlob=Patch_Boundary_Edge{nn};
    %E_InternalGlob=setdiff(E_Glob,E_BoundaryGlob);
    E_InternalGlob=E_Glob;
    
    
    
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
%     b_sigma1Loc= - A(E_InternalGlob1,     Tot_ExternalGlob) * x(Tot_ExternalGlob);
%     b_sigma2Loc= - A(E_InternalGlob2,     Tot_ExternalGlob) * x(Tot_ExternalGlob);
%     b_disp1Loc=  - A(N_InternalGlob1,     Tot_ExternalGlob) * x(Tot_ExternalGlob);
%     b_disp2Loc=  - A(N_InternalGlob2,     Tot_ExternalGlob) * x(Tot_ExternalGlob);    
%     b_sigma1Loc=b_sigma1Loc+b_sigma1(E_InternalGlob);
%     b_sigma2Loc=b_sigma2Loc+b_sigma2(E_InternalGlob);
%     b_disp1Loc=b_disp1Loc+b_disp1(N_InternalGlob);
%     b_disp2Loc=b_disp2Loc+b_disp2(N_InternalGlob);    
    b_sigma1Loc= b_sigma1(E_InternalGlob)- A(E_InternalGlob1,     :) *x ;
    b_sigma2Loc= b_sigma2(E_InternalGlob)- A(E_InternalGlob2,     :) *x ;
    b_disp1Loc= b_disp1(N_InternalGlob)  - A(N_InternalGlob1,     :) *x ;
    b_disp2Loc= b_disp2(N_InternalGlob)  - A(N_InternalGlob2,     :) *x ;    
  
    
%     A_Loc{3,1}(N_BoundaryLoc,:)=0; 
%     A_Loc{3,2}(N_BoundaryLoc,:)=0; 
%     A_Loc{3,3}(N_BoundaryLoc,:)=0;
%     A_Loc{3,4}(N_BoundaryLoc,:)=0;      
%     A_Loc{4,1}(N_BoundaryLoc,:)=0; 
%     A_Loc{4,2}(N_BoundaryLoc,:)=0; 
%     A_Loc{4,3}(N_BoundaryLoc,:)=0; 
%     A_Loc{4,4}(N_BoundaryLoc,:)=0; 
    
%     cont=0;
%     for nn_loc=N_BoundaryLoc
%         cont=cont+1;
%         A_Loc{3,3}(nn_loc,nn_loc)=1;
%         A_Loc{4,4}(nn_loc,nn_loc)=1;
%         b_disp1Loc(nn_loc)=x(2*NE+N_BoundaryGlob(cont));
%         b_disp2Loc(nn_loc)=x(2*NE+N+N_BoundaryGlob(cont));
%     end


    % since we express external force just through the displacement field
    % on the boundary we put zero forces (IS THIS CORRECT OR UOT?)

 
    
    
    b_Loc=[b_sigma1Loc;b_sigma2Loc;b_disp1Loc;b_disp2Loc];
   
    
    
    
    
    Tot_otherGlob=setdiff(1:length(b),Tot_ExternalGlob); 
    
    
    
   
% by construction the node nn is on GammaC
% N_contact=intersect(N_Glob,mesh.N_contact);
% E_contact=intersect(mesh.E_bc_contact,E_InternalGlob1);
N_contact=intersect(N_InternalGlob,mesh.N_contact);
E_contact=intersect(mesh.E_bc_contact,E_InternalGlob1);
[jj,nn]

E_contact_Loc=[];
N_contact_Loc=[];
for kk=1:length(E_contact)
    E_contact_Loc(kk,1)=mesh.E_contact_map(E_contact(kk));
end



E_loc_map=containers.Map(E_InternalGlob,1:length(E_InternalGlob));
mesh_Loc.E_contact=[];
for kk=1:length(E_contact)
    mesh_Loc.E_contact(kk)=E_loc_map(E_contact(kk));
end

N_loc_map=containers.Map(N_contact,1:length(N_contact));
mesh_Loc.N_contact=[];


for kk=1:length(N_contact)
    mesh_Loc.N_contact(kk)=N_loc_map(N_contact(kk));
    N_contact_Loc(kk)=mesh.N_contact_map(N_contact(kk));
end


Tot_InternalGlobContact=[E_InternalGlob,E_InternalGlob+NE,N_contact+2*NE,N_contact+2*NE+N];
actual_solution=x(Tot_InternalGlobContact);    
    
    
    
    
mesh_Loc.NE=length(E_InternalGlob); 
mesh_Loc.N=length(N_InternalGlob);


Constraint_Loc.CheckConstraintE1=Constraint.CheckConstraintE1(E_contact_Loc,Tot_InternalGlob);
Constraint_Loc.RhsE1=Constraint.RhsE1(E_contact_Loc,1) - Constraint.CheckConstraintE1(E_contact_Loc,:)*x;

Constraint_Loc.CheckConstraintN1=Constraint.CheckConstraintN1(N_contact_Loc,Tot_InternalGlob);
Constraint_Loc.RhsN1=Constraint.RhsN1(N_contact_Loc,1) - Constraint.CheckConstraintN1(N_contact_Loc,:)*x;

Constraint_Loc.CheckConstraintE2=Constraint.CheckConstraintE2(E_contact_Loc,Tot_InternalGlob);
Constraint_Loc.RhsE2=Constraint.RhsE2(E_contact_Loc,1) - Constraint.CheckConstraintE2(E_contact_Loc,:)*x;

Constraint_Loc.CheckConstraintN2=Constraint.CheckConstraintN2(N_contact_Loc,Tot_InternalGlob);
Constraint_Loc.RhsN2=Constraint.RhsN2(N_contact_Loc,1) - Constraint.CheckConstraintN2(N_contact_Loc,:)*x;


% Constraint_Loc.CheckConstraintE1=Constraint.CheckConstraintE1(E_contact_Loc,Tot_InternalGlob);
% Constraint_Loc.RhsE1=Constraint.RhsE1(E_contact_Loc,1) - Constraint.CheckConstraintE1(E_contact_Loc,Tot_ExternalGlob)*x(Tot_ExternalGlob);
% 
% Constraint_Loc.CheckConstraintN1=Constraint.CheckConstraintN1(N_contact_Loc,Tot_InternalGlob);
% Constraint_Loc.RhsN1=Constraint.RhsN1(N_contact_Loc,1) - Constraint.CheckConstraintN1(N_contact_Loc,Tot_ExternalGlob)*x(Tot_ExternalGlob);
% 
% Constraint_Loc.CheckConstraintE2=Constraint.CheckConstraintE2(E_contact_Loc,Tot_InternalGlob);
% Constraint_Loc.RhsE2=Constraint.RhsE2(E_contact_Loc,1) - Constraint.CheckConstraintE2(E_contact_Loc,Tot_ExternalGlob)*x(Tot_ExternalGlob);
% 
% Constraint_Loc.CheckConstraintN2=Constraint.CheckConstraintN2(N_contact_Loc,Tot_InternalGlob);
% Constraint_Loc.RhsN2=Constraint.RhsN2(N_contact_Loc,1) - Constraint.CheckConstraintN2(N_contact_Loc,Tot_ExternalGlob)*x(Tot_ExternalGlob);
% 

%intersect(E_InternalGlob,Constraint.WorkingSetE)

Constraint_Loc.WorkingSetE=[];
tmpE=intersect(E_InternalGlob,Constraint.WorkingSetE);
for pp=1:length(tmpE)
Constraint_Loc.WorkingSetE(pp)=E_loc_map(tmpE(pp));
end

Constraint_Loc.WorkingSetN=[];
tmpN=intersect(N_InternalGlob,Constraint.WorkingSetN);
for pp=1:length(tmpN)
    %if lengt(tmpN)>0, then
Constraint_Loc.WorkingSetN(pp)=1;
end

% now we redefine
%mesh_Loc.N_contact=1;

is_the_right_solution=false;

    while(is_the_right_solution==false)
    [c_Loc,is_the_right_solution,Constraint_Loc] = activeset3(A_Loc,b_Loc,Constraint_Loc,mesh_Loc);
    end
    
    % we remove the values of local E_contact, N_contact
    Constraint.WorkingSetE=setdiff(Constraint.WorkingSetE,tmpE);
    Constraint.WorkingSetN=setdiff(Constraint.WorkingSetN,tmpN);    
    % then, after the local activeset, we add the new ones
    tmpWorkingSetE=[Constraint.WorkingSetE,E_InternalGlob(Constraint_Loc.WorkingSetE)];
    Constraint.WorkingSetE=unique(tmpWorkingSetE);
    tmpWorkingSetN=[Constraint.WorkingSetN,N_InternalGlob(Constraint_Loc.WorkingSetN)];
    Constraint.WorkingSetN=unique(tmpWorkingSetN);
    
    x(Tot_InternalGlob')=x(Tot_InternalGlob')+c_Loc;
    
   
%     sol=M_Normal_Tangent*x;
%     print_displacement_solution(meshprova,sol(1+2*meshprova{L}.NE:2*meshprova{L}.NE+meshprova{L}.N)',sol(1+2*meshprova{L}.NE+meshprova{L}.N:end)');

end





end


end

