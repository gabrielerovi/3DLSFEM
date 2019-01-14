function [A,b,A11,A12,A13,A14,A22,A23,A24,A33,A34,A44]= A_with_BC(A, b,A11,A12,A13,A14,A22,A23,A24,A33,A34,A44,mesh,contact)

E_label=mesh.E_label;
N_label=mesh.N_label;

N=mesh.N;
NE=mesh.NE;
E_remove1=mesh.E_remove;
E_remove=E_remove1;
E_remove2=E_remove1+NE;
N_remove1=mesh.N_remove;
N_remove=N_remove1;
N_remove1=N_remove1+2*NE;
N_remove2=N_remove1+N;
Remove=[E_remove1 E_remove2 N_remove1 N_remove2];


sigma1=b(1:NE);
sigma2=b(NE+1:2*NE);
disp1=b(2*NE+1:2*NE+N);
disp2=b(2*NE+N+1:end);

A11(E_remove,:)=0; A12(E_remove,:)=0; A13(E_remove,:)=0; A14(E_remove,:)=0;
A22(E_remove,:)=0; A23(E_remove,:)=0; A24(E_remove,:)=0;
A33(N_remove,:)=0; A34(N_remove,:)=0; A44(N_remove,:)=0;
EcontBC=0;
type_of_dof=2;
for ii=E_remove
    U_xy=[0;0];
    EcontBC=EcontBC+1;
    tmp=add_boundary_bc_elasticity2D(U_xy,type_of_dof,E_label(EcontBC), contact);  
    coeff = RT_dirichlet_coeff(E_remove(EcontBC), mesh);
    sigma1(ii)=tmp(1)/coeff;
    sigma2(ii)=tmp(2)/coeff;
    A11(ii,ii)=1;  
    A22(ii,ii)=1; 
end

NcontBC=0;
type_of_dof=1;
for ii=N_remove
    U_xy=[0;0];
    NcontBC=NcontBC+1;
    tmp=add_boundary_bc_elasticity2D(U_xy,type_of_dof,N_label(NcontBC), contact);   
    disp1(ii)=tmp(1);
    disp2(ii)=tmp(2);
    A33(ii,ii)=1;  
    A44(ii,ii)=1; 

end




A(Remove,:)=0;
for ii=Remove
   A(ii,ii)=1;
end







contact=0;
EcontBC=0;
type_of_dof=2;


contact=0;
NcontBC=0;   
type_of_dof=1;


b=[sigma1;sigma2;disp1;disp2];

end