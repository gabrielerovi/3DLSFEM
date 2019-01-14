function U=add_boundary_bc_elasticity3D(type_of_dof,label,dim)

[dirichlet]= boundary_value_bool(type_of_dof);
U=dirichlet(label,1:dim)';


end