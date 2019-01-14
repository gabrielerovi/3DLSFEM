function U=add_boundary_bc_elasticity2D(U_xy,type_of_dof,label,contact)



% bool is a vector that tells us if the displacement is constrained:
% only in the imposed-displacement direction -> roller in the orthogonal one
% totally -> hinge with imposed displacement
[dirichlet,n_and_or_t,bool_bc]= boundary_value_bool(type_of_dof);
U_bc=dirichlet(label,[1 2]);
U_bc=U_bc';

if(contact==0)
U=U_bc;
else
bool_n_t=n_and_or_t(label,[1,2]);
bool_n_t=bool_n_t';
% define the computed vector in the coordinate system of U_bc
% now the first basis vector= U_bc 
% this means that only this component
unit=U_bc/norm(U_bc);
Q=[ unit(1)  unit(2);
   -unit(2)  unit(1)];

U_nt=Q*U_xy;
U_bc=Q*U_bc;
U=Q'*(U_nt.*(1-bool_n_t)+ U_bc .* bool_n_t)
end
end