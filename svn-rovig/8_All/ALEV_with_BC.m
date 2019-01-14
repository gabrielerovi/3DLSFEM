function [Ass,Asu,Aus,Auu,bs,bu,sigma1,sigma2,disp1,disp2,...
            A11,A12,A13,A14,...
            A21,A22,A23,A24,...
            A31,A32,A33,A34,...
            A41,A42,A43,A44]= A_with_BC(b,A11,A12,A13,A14,A22,A23,A24,A33,A34,A44,grid,contact)


L=size(grid);
L=L(1);
    

lev=L;
mesh=grid{lev};
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

A21{lev}=A12{lev}';
A31{lev}=A13{lev}';
A32{lev}=A23{lev}';
A41{lev}=A14{lev}';
A42{lev}=A24{lev}';
A43{lev}=A34{lev}';




A11{lev}(E_remove,:)=0; A12{lev}(E_remove,:)=0; A13{lev}(E_remove,:)=0; A14{lev}(E_remove,:)=0;
A21{lev}(E_remove,:)=0; A22{lev}(E_remove,:)=0; A23{lev}(E_remove,:)=0; A24{lev}(E_remove,:)=0;
A31{lev}(N_remove,:)=0; A32{lev}(N_remove,:)=0; A33{lev}(N_remove,:)=0; A34{lev}(N_remove,:)=0; 
A41{lev}(N_remove,:)=0; A42{lev}(N_remove,:)=0; A43{lev}(N_remove,:)=0; A44{lev}(N_remove,:)=0;


EcontBC=0;
type_of_dof=2;
for ii=E_remove
    U_xy=[0;0];
    EcontBC=EcontBC+1;
    tmp=add_boundary_bc_elasticity2D(U_xy,type_of_dof,E_label(EcontBC), contact);  
    coeff = RT_dirichlet_coeff(E_remove(EcontBC), mesh);
    sigma1(ii)=tmp(1)/coeff;
    sigma2(ii)=tmp(2)/coeff;
    A11{lev}(ii,ii)=1;  
    A22{lev}(ii,ii)=1; 
end

NcontBC=0;
type_of_dof=1;
for ii=N_remove
    U_xy=[0;0];
    NcontBC=NcontBC+1;
    tmp=add_boundary_bc_elasticity2D(U_xy,type_of_dof,N_label(NcontBC), contact);   
    disp1(ii)=tmp(1);
    disp2(ii)=tmp(2);
    A33{lev}(ii,ii)=1;  
    A44{lev}(ii,ii)=1; 

end

Ass{lev}=[A11{lev} A12{lev};
          A21{lev} A22{lev}];
Asu{lev}=[A13{lev} A14{lev};
          A23{lev} A24{lev}];
Aus{lev}=[A31{lev} A32{lev};
          A41{lev} A42{lev}];      
Auu{lev}=[A33{lev} A34{lev};
          A43{lev} A44{lev}];
      
      
bs=[sigma1;sigma2];
bu=[disp1;disp2];





















for lev=1:L-1
    
A21{lev}=A12{lev}';
A31{lev}=A13{lev}';
A32{lev}=A23{lev}';
A41{lev}=A14{lev}';
A42{lev}=A24{lev}';
A43{lev}=A34{lev}';

mesh=grid{lev};

N=mesh.N;
NE=mesh.NE;
E_remove=mesh.E_remove;
N_remove=mesh.N_remove;

A11{lev}(E_remove,:)=0; A12{lev}(E_remove,:)=0; A13{lev}(E_remove,:)=0; A14{lev}(E_remove,:)=0;
A21{lev}(E_remove,:)=0; A22{lev}(E_remove,:)=0; A23{lev}(E_remove,:)=0; A24{lev}(E_remove,:)=0;
A31{lev}(N_remove,:)=0; A32{lev}(N_remove,:)=0; A33{lev}(N_remove,:)=0; A34{lev}(N_remove,:)=0; 
A41{lev}(N_remove,:)=0; A42{lev}(N_remove,:)=0; A43{lev}(N_remove,:)=0; A44{lev}(N_remove,:)=0;

for ii=E_remove
    A11{lev}(ii,ii)=1;  
    A22{lev}(ii,ii)=1; 
end

for ii=N_remove 
    A33{lev}(ii,ii)=1;  
    A44{lev}(ii,ii)=1; 

end


Ass{lev}=[A11{lev} A12{lev};
          A21{lev} A22{lev}];
Asu{lev}=[A13{lev} A14{lev};
          A23{lev} A24{lev}];
Aus{lev}=[A31{lev} A32{lev};
          A41{lev} A42{lev}];      
Auu{lev}=[A33{lev} A34{lev};
          A43{lev} A44{lev}];
             
end









end