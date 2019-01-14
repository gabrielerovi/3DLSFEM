function [mesh,h]=mesh(parameters)

COARSE=parameters.COARSE;
number_of_levels=parameters.number_of_levels;
dim=parameters.dim;
N_components=parameters.N_components;
E_components=parameters.E_components;
F_components=parameters.F_components;


% generate L nested meshes
mesh =create_meshes(COARSE,dim,N_components,E_components,F_components,parameters);
mesh =create_meshes_coarser_to_finer(mesh,number_of_levels,dim);
mesh=node_and_face_normal_boundary(mesh,parameters);
mesh=mesh_contact_dofs(mesh,parameters);
h=meshwidth(mesh);

end