function  mesh =creaste_meshes_coarser_to_fines(old_mesh,L,dim)
mesh=cell(L,1);
old_L=length(old_mesh);


% mesh{1} is equal to
mesh{1}=old_mesh{old_L};
% mesh{1}.node_per_elem=old_mesh{old_L}.node_per_elem;
% mesh{1}.edge_per_elem=old_mesh{old_L}.edge_per_elem;
% mesh{1}.face_per_elem=old_mesh{old_L}.face_per_elem;
% boundary=old_mesh{old_L}.boundary;
% mesh{1}.node_per_elem=old_mesh{old_L}.node_per_elem;
% mesh{1}.edge_per_elem=old_mesh{old_L}.edge_per_elem;
% mesh{1}.face_per_elem=old_mesh{old_L}.face_per_elem;
% mesh{1}.node=old_mesh{old_L}.node;
% % mesh{1}.elem=sort(old_mesh{old_L}.elem;,2);
% mesh{1}.boundary=boundary;
% mesh{1}.N=length(node(:,1));
% mesh{1}.NT=length(elem(:,1));
% mesh{1}.N_components=N_components;
% mesh{1}.E_components=E_components;
% mesh{1}.F_components=F_components;


% N=length(node(:,1));
% NT=length(elem(:,1));


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
 [mesh{1}.NE,mesh{1}.edge,mesh{1}.elemE,mesh{1}.E_to_T,mesh{1}.N_to_T]=create_edge_structures3D(mesh{1}.edge_per_elem,mesh{1}.node_per_elem,mesh{1}.node,mesh{1}.elem);
 [mesh{1}.NF,mesh{1}.face,mesh{1}.elemF,mesh{1}.F_to_T]=create_face_structures3D(mesh{1}.face_per_elem,mesh{1}.node,mesh{1}.elem);
% this function is called only for the roguh mesh
% then for each other fine mesh L, T_face_bc(L) is built based on T_face_bc(L-1)
mesh{1}.T_face_bc=T_to_bc3D(mesh{1}.elem,mesh{1}.boundary,mesh{1}.face_per_elem);

 [mesh{1}.flag_for_bc_face,mesh{1}.flag_for_each_face] = flag_boundary_faces3D  (mesh{1}.boundary,mesh{1}.face);
[mesh{1}.N_bc,mesh{1}.F_bc,mesh{1}.T_bc]=boundary_flags3D(mesh{1}.boundary, mesh{1}.N,  mesh{1}.face,mesh{1}.F_to_T,mesh{1}.face_per_elem);  

[mesh{1}.N_dirichlet,mesh{1}.N_label ]=is_surface_dirichlet(mesh{1}.N_bc,1); 
% [mesh{1}.E_dirichlet,mesh{1}.E_label]=is_surface_dirichlet(mesh{1}.E_bc,2); 
[mesh{1}.F_dirichlet,mesh{1}.F_label]=is_surface_dirichlet(mesh{1}.F_bc,3); 
mesh{1}.N_remove=find(mesh{1}.N_dirichlet); 
% mesh{1}.E_remove=find(mesh{1}.E_dirichlet); 
mesh{1}.F_remove=find(mesh{1}.F_dirichlet); 

for lev=2:L
mesh{lev}.N_components=mesh{lev-1}.N_components;
mesh{lev}.E_components=mesh{lev-1}.E_components;
mesh{lev}.F_components=mesh{lev-1}.F_components;
mesh{lev}.node_per_elem=mesh{lev-1}.node_per_elem;
mesh{lev}.edge_per_elem=mesh{lev-1}.edge_per_elem;
mesh{lev}.face_per_elem=mesh{lev-1}.face_per_elem;

[mesh{lev}.node,mesh{lev}.elem,mesh{lev}.T_face_bc,mesh{lev}.boundary] = uniformrefine3(mesh{lev-1}.node,mesh{lev-1}.elem,mesh{lev-1}.T_face_bc);
mesh{lev}.NT=length(mesh{lev}.elem(:,1));
mesh{lev}.N=length(mesh{lev}.node(:,1));

 [mesh{lev}.NE,mesh{lev}.edge,mesh{lev}.elemE,mesh{lev}.E_to_T,mesh{lev}.N_to_T]=... 
create_edge_structures3D(mesh{lev}.edge_per_elem,mesh{lev}.node_per_elem,mesh{lev}.node,mesh{lev}.elem);
 [mesh{lev}.NF,mesh{lev}.face,mesh{lev}.elemF,mesh{lev}.F_to_T]= ...
create_face_structures3D(mesh{lev}.face_per_elem,mesh{lev}.node,mesh{lev}.elem);
 mesh{lev}.N=length(mesh{lev}.node(:,1));
 [mesh{lev}.flag_for_bc_face,mesh{lev}.flag_for_each_face] = ... 
flag_boundary_faces3D  (mesh{lev-1}.boundary,mesh{lev-1}.face);
 [mesh{lev}.N_bc,mesh{lev}.F_bc,mesh{lev}.T_bc]=boundary_flags3D(mesh{lev}.boundary, mesh{lev}.N, mesh{lev}.face,mesh{lev}.F_to_T,mesh{lev}.face_per_elem);    

[mesh{lev}.N_dirichlet,mesh{lev}.N_label]=is_surface_dirichlet(mesh{lev}.N_bc,1); 
% [mesh{lev}.E_dirichlet,mesh{lev}.E_label]=is_surface_dirichlet(mesh{lev}.E_bc,2); 
[mesh{lev}.F_dirichlet,mesh{lev}.F_label]=is_surface_dirichlet(mesh{lev}.F_bc,3); 
mesh{lev}.N_remove=find(mesh{lev}.N_dirichlet); 
% mesh{lev}.E_remove=find(mesh{lev}.E_dirichlet); 
mesh{lev}.F_remove=find(mesh{lev}.F_dirichlet); 
end

end


end