function  mesh =create_meshes(L,dim,N_components,E_components,F_components,parameters)
mesh=cell(L,1);
if(dim==2)
node_per_elem=3;
face_per_elem=3;
edge_per_elem=3;
% 2D example
if(parameters.create_only_square==true)
a=parameters.a;
b=parameters.b;
node = [0,0; a,0; a,b; 0,b];
elem = [2,3,1; 4,1,3];
boundary=[1 2 1; 2 3 2; 3 4 3; 4 1 4];
else
[node,elem,boundary]=mesh_read(parameters.input_mesh);
elem=elem(:,[1,2,3]);
end
% node = [0,0; a,0; a,b; 0,b];
% elem =[    2,3,1; 4,1,3];
% boundary=[1 2 1; 2 3 2; 3 4 3; 4 1 4];




%  node = [0.0 0.0; 0.5 0.0; 1.0 0.0;
%          0.0 0.5; 0.5 0.5; 1.0 0.5;
%          0.0 1.0; 0.5 1.0; 1.0 1.0;]
%  elem = [
% 1 2 5 ;
% 1 5 4 ;
% 2 3 6 ;
% 2 6 5 ;
% 4 5 8 ;
% 4 8 7 ;
% 5 6 9 ;
% 5 9 8 ;
% ];
%  boundary=[1 2 1; 2 3 1; 
%            3 6 2; 6 9 2; 
%            9 8 3; 8 7 3; 
%            7 4 4; 4 1 4];









% 
% a=0.5;
% b=0.2;
% node = [0,0; a,0; a,b; 0,b; 2*a,0;2*a,b];
% elem =[    2,3,1; 4,1,3; 2 5 6; 2 6 3;];
% boundary=[1 2 1; 2 5 1; 5 6 2; 6 3 3; 3 4 3; 1 4 4;];


% node = [0,0; 1,0; 1,1; 0,1];
% node = [ 0,  0; a/4,0; a/2,0; 3/4*a,0; a,0; 
%          0,  b; a/4,b; a/2,b; 3/4*a,b; a,b;];
% elem = [1,6,7;1,2,7;2,7,8;2,3,8;3,8,9;3,4,9;4,10,9;4,5,10];
% boundary=[1 2 1; 2 3 1; 3 4 1; 4 5 1; ...
%           5 10 2; ...
%           6 7 3; 7 8 3; 8 9 3; 9 10 3; ...
%           6 1 4;];


% 8     9     10
% 4   5   6   7
% 1     2     3
% node = [0,0; a/2,0; a,0; 
%         0,b/2; a/3,b/2; 2/3*a,b/2; a,b/2;
%         0,b; a/2,b; a,b;];
% elem=[1 2 5; 2 5 6;   2 6 3; 1 4 5;  3 6 7;
%       8 9 5; 9 10 6; 5 6 9;  8 4 5; 10 6 7;]
% 
% boundary=[ 1 2 1; 2 3 1; 3 7 2; 7 10 2; 10 9 3; 8 9 3; 1 4 4; 4 8 4;];




% 
%   node = [0,0; a,0; a,b; 0,b;a/2,b/2];
% elem = [1,2,5; 2,3,5;3,4,5;4,1,5];
% boundary=[1 2 1; 2 3 2; 3 4 3; 4 1 4];



% 
% 
%  node = [       0         0;
%     1.0000         0;
%     1.0000    1.0000;
%          0    1.0000;
%     0.5000         0;
%     0.5000    0.5000;
%          0    0.5000;
%     1.0000    0.5000;
%     0.5000    1.0000];
%  elem = [
%      1     5     6;
%      2     5     8;
%      3     6     8;
%      5     6     8;
%      1     6     7;
%      3     6     9;
%      4     7     9;
%      6     7     9;
% ];
%  boundary=[1 5 1; 5 2 1; 
%            2 8 2; 8 3 2; 
%            3 9 3; 9 4 3; 
%            4 7 4; 7 1 4];

% 
% node = [0,0; 1,0; 1,1; 0,1; 1,2; 0,2; 2,1;2,2];
% elem = [2,3,1; 4,1,3; 3 4 5; 4 5 6; 3 7 8; 3 5 8];
% boundary=[1 2 1; 2 3 2; 3 7 3; 7 8 4; 8 5 5; 5 6 5; 6 4 6; 4 1 4];
else
node_per_elem=4;
face_per_elem=4;
edge_per_elem=6;
% 3D example
node = [-1,-1,-1; 1,-1,-1; 1,1,-1; -1,1,-1; -1,-1,1; 1,-1,1; 1,1,1; -1,1,1];
elem = [1 2 3 7; 1 4 3 7; 1 5 6 7; 1 5 8 7; 1 2 6 7; 1 4 8 7];
boundary=[1 2  3 100; 2 3 7 2100; 1 2 7 3100; 1 3 7 4100];
node = [0,0,0; 1,0,0; 0,1,0; 0,0,1; ];
elem = [1 2 3 4; ];
boundary=[1 2  3 10; 1 2 4 20; 1 3 4 30; 2 3 4 40];
end
mesh{1}.node_per_elem=node_per_elem;
mesh{1}.edge_per_elem=edge_per_elem;
mesh{1}.face_per_elem=face_per_elem;
mesh{1}.node=node;
mesh{1}.elem=sort(elem,2);
mesh{1}.boundary=boundary;
mesh{1}.N=length(node(:,1));
mesh{1}.NT=length(elem(:,1));
mesh{1}.N_components=N_components;
mesh{1}.E_components=E_components;
mesh{1}.F_components=F_components;


N=length(node(:,1));
NT=length(elem(:,1));


if(dim==2)
[mesh{1}.NE,mesh{1}.edge,mesh{1}.elemE,mesh{1}.E_to_T,mesh{1}.N_to_T]= ... 
 create_edge_structures2D(mesh{1}.edge_per_elem,mesh{1}.node,mesh{1}.elem);
mesh{1}.T_edge_bc=T_to_bc2D(mesh{1}.elem,mesh{1}.boundary,mesh{1}.edge_per_elem);
[mesh{1}.flag_for_bc_edge,mesh{1}.flag_for_each_edge] = ... 
flag_boundary_edges2D  (mesh{1}.boundary,mesh{1}.edge);
[mesh{1}.N_bc,mesh{1}.E_bc,mesh{1}.T_bc]=boundary_flags2D ...
(mesh{1}.boundary, mesh{1}.N, mesh{1}.edge,mesh{1}.E_to_T);

[mesh{1}.N_dirichlet,mesh{1}.N_label] =is_surface_dirichlet(mesh{1}.N_bc,1); 
[mesh{1}.E_dirichlet,mesh{1}.E_label] =is_surface_dirichlet(mesh{1}.E_bc,2); 
mesh{1}.N_remove=find(mesh{1}.N_dirichlet); 
mesh{1}.E_remove=find(mesh{1}.E_dirichlet); 
[mesh{1}.ND_nodes_on_edge_boundary,mesh{1}.ND_nodes_reorder,mesh{1}.ND_lengths]=nodes_on_edge_boundary(mesh{1});

for lev=2:L
mesh{lev}.node_per_elem=node_per_elem;
mesh{lev}.edge_per_elem=edge_per_elem;
mesh{lev}.face_per_elem=face_per_elem;


mesh{lev}.edge_per_elem=edge_per_elem;
mesh{lev}.face_per_elem=face_per_elem;

[mesh{lev}.node,mesh{lev}.elem,mesh{lev}.T_edge_bc,mesh{lev}.boundary] = ...
    uniformrefine(mesh{lev-1}.node,mesh{lev-1}.elem,mesh{lev-1}.T_edge_bc);

mesh{lev}.edge_per_elem=edge_per_elem;
mesh{lev}.face_per_elem=face_per_elem;
mesh{lev}.N=length(mesh{lev}.node(:,1));
mesh{lev}.NT=length(mesh{lev}.elem(:,1));

[mesh{lev}.NE,mesh{lev}.edge,mesh{lev}.elemE,mesh{lev}.E_to_T,mesh{lev}.N_to_T]= ... 
    create_edge_structures2D(mesh{lev}.edge_per_elem,mesh{lev}.node,mesh{lev}.elem);
[mesh{lev}.flag_for_bc_edge,mesh{lev}.flag_for_each_edge] = ...
    flag_boundary_edges2D  (mesh{lev}.boundary,mesh{lev}.edge);

 [mesh{lev}.N_bc,mesh{lev}.E_bc,mesh{lev}.T_bc]=boundary_flags2D ... 
     (mesh{lev}.boundary, mesh{lev}.N, mesh{lev}.edge,mesh{lev}.E_to_T); 

[mesh{lev}.N_dirichlet,mesh{lev}.N_label]=is_surface_dirichlet(mesh{lev}.N_bc,1); 
[mesh{lev}.E_dirichlet,mesh{lev}.E_label]=is_surface_dirichlet(mesh{lev}.E_bc,2); 
mesh{lev}.N_remove=find(mesh{lev}.N_dirichlet); 
mesh{lev}.E_remove=find(mesh{lev}.E_dirichlet); 

[mesh{lev}.ND_nodes_on_edge_boundary,mesh{lev}.ND_nodes_reorder,mesh{lev}.ND_lengths]=nodes_on_edge_boundary(mesh{lev});


end


else
 [mesh{1}.NE,mesh{1}.edge,mesh{1}.elemE,mesh{1}.E_to_T]=create_edge_structures3D(mesh{1}.edge_per_elem,mesh{1}.node,mesh{1}.elem);
 [mesh{1}.NF,mesh{1}.face,mesh{1}.elemF,mesh{1}.F_to_T]=create_face_structures3D(mesh{1}.face_per_elem,mesh{1}.node,mesh{1}.elem);
% this function is called only for the roguh mesh
% then for each other fine mesh L, T_face_bc(L) is built based on T_face_bc(L-1)
mesh{1}.T_face_bc=T_to_bc3D(mesh{1}.elem,mesh{1}.boundary,mesh{1}.face_per_elem);

 [mesh{1}.flag_for_bc_face,mesh{1}.flag_for_each_face] = flag_boundary_faces3D  (mesh{1}.boundary,mesh{1}.face);
 %[mesh{1}.N_bc,mesh{1}.E_bc,mesh{1}.F_bc,mesh{1}.T_bc]=boundary_flags3D(mesh{1}.boundary, mesh{1}.N, mesh{1}.edge, mesh{1}.face,mesh{1}.F_to_T);  

[mesh{1}.N_dirichlet,mesh{1}.N_label ]=is_surface_dirichlet(mesh{1}.N_bc,1); 
[mesh{1}.E_dirichlet,mesh{1}.E_label]=is_surface_dirichlet(mesh{1}.E_bc,2); 
[mesh{1}.F_dirichlet,mesh{1}.F_label]=is_surface_dirichlet(mesh{1}.F_bc,3); 
mesh{1}.N_remove=find(mesh{1}.N_dirichlet); 
mesh{1}.E_remove=find(mesh{1}.E_dirichlet); 
mesh{1}.F_remove=find(mesh{1}.F_dirichlet); 

for lev=2:L
mesh{lev}.node_per_elem=node_per_elem;
mesh{lev}.edge_per_elem=edge_per_elem;
mesh{lev}.face_per_elem=face_per_elem;

[mesh{lev}.node,mesh{lev}.elem,mesh{lev}.T_face_bc,mesh{lev}.boundary] = uniformrefine3(mesh{lev-1}.node,mesh{lev-1}.elem,mesh{lev-1}.T_face_bc);

 [mesh{lev}.NE,mesh{lev}.edge,mesh{lev}.elemE,mesh{lev}.E_to_T]=... 
create_edge_structures3D(mesh{lev}.edge_per_elem,mesh{lev}.node,mesh{lev}.elem);
 [mesh{lev}.NF,mesh{lev}.face,mesh{lev}.elemF,mesh{lev}.F_to_T]= ...
create_face_structures3D(mesh{lev}.face_per_elem,mesh{lev}.node,mesh{lev}.elem);
 mesh{lev}.N=length(mesh{lev}.node(:,1));
 [mesh{lev}.flag_for_bc_face,mesh{lev}.flag_for_each_face] = ... 
flag_boundary_faces3D  (mesh{lev-1}.boundary,mesh{lev-1}.face);
%  [mesh{lev}.N_bc,mesh{lev}.E_bc,mesh{lev}.F_bc,mesh{lev}.T_bc]= ...
% boundary_flags3D(mesh{lev-1}.boundary, mesh{lev-1}.N, mesh{lev-1}.edge, mesh{lev-1}.face,mesh{lev-1}.F_to_T);  

[mesh{lev}.N_dirichlet,mesh{lev}.N_label]=is_surface_dirichlet(mesh{lev}.N_bc,1); 
[mesh{lev}.E_dirichlet,mesh{lev}.E_label]=is_surface_dirichlet(mesh{lev}.E_bc,2); 
[mesh{lev}.F_dirichlet,mesh{lev}.F_label]=is_surface_dirichlet(mesh{lev}.F_bc,3); 
mesh{lev}.N_remove=find(mesh{lev}.N_dirichlet); 
mesh{lev}.E_remove=find(mesh{lev}.E_dirichlet); 
mesh{lev}.F_remove=find(mesh{lev}.F_dirichlet); 
end

end

end