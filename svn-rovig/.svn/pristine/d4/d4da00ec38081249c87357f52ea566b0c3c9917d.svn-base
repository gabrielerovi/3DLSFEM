function [WorkingSetNormal_E,WorkingSetNormal_N]=WorkingSetsNormalTangent(x,mesh)

toll=10^(-12);

E_contact=mesh.E_contact;
N_contact=mesh.N_contact;

NE=mesh.NE;
N=mesh.N;

% sigma_normal=sparse(NE,1);
% sigma_tangent=sparse(NE,1);
% disp_normal=sparse(N,1);
% disp_tangent=sparse(N,1);
% g_normal=sparse(N,1)

sigma_normal  = x(E_contact);
sigma_tangent = x(E_contact + NE );
disp_normal   = x(N_contact + 2*NE );
disp_tangent  = x(N_contact + N + 2 * NE );


g_normal=gap_function(mesh.node(N_contact,1),mesh.node(N_contact,2));

WorkingSetNormal_E=E_contact(find(abs(sigma_normal)<toll));
WorkingSetNormal_N=N_contact(find(abs(disp_normal -g_normal)<toll));



end