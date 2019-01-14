function [Ass,b,A11_lev,A12_lev,A21_lev,A22_lev ]=A_meshwidth_with_BC(RTCtoRTF,b,qrule,alpha,beta,coeff_equilibrium,coeff_constitutive,coeff_symmetry,mesh,contact)


L=size(mesh);
L=L(1);
    
h=meshwidth(mesh);

%h=ones(L,1);
lev=L;
A11Eq_lev{lev}=assembling2DRTRT(mesh{L},qrule,alpha,beta,1,1,coeff_equilibrium,0,0);
A12Eq_lev{lev}=assembling2DRTRT(mesh{L},qrule,alpha,beta,1,2,coeff_equilibrium,0,0);
A21Eq_lev{lev}=A12Eq_lev{lev}';
A22Eq_lev{lev}=assembling2DRTRT(mesh{L},qrule,alpha,beta,2,2,coeff_equilibrium,0,0);

A11AsymConst_lev{lev}=assembling2DRTRT(mesh{L},qrule,alpha,beta,1,1,0,coeff_constitutive,coeff_symmetry);
A12AsymConst_lev{lev}=assembling2DRTRT(mesh{L},qrule,alpha,beta,1,2,0,coeff_constitutive,coeff_symmetry);
A21AsymConst_lev{lev}=A12AsymConst_lev{lev}';
A22AsymConst_lev{lev}=assembling2DRTRT(mesh{L},qrule,alpha,beta,2,2,0,coeff_constitutive,coeff_symmetry);



if(L>1)
for lev=L-1:-1:1   
    
    A11Eq_lev{lev}=h(lev)^2.* RTCtoRTF{lev}'*A11Eq_lev{lev+1}*RTCtoRTF{lev};
    A12Eq_lev{lev}=h(lev)^2.* RTCtoRTF{lev}'*A12Eq_lev{lev+1}*RTCtoRTF{lev};
    A21Eq_lev{lev}=A12Eq_lev{lev}';
    A22Eq_lev{lev}=h(lev)^2.* RTCtoRTF{lev}'*A22Eq_lev{lev+1}*RTCtoRTF{lev};
 
    A11AsymConst_lev{lev}=RTCtoRTF{lev}'*A11AsymConst_lev{lev+1}*RTCtoRTF{lev};
    A12AsymConst_lev{lev}=RTCtoRTF{lev}'*A12AsymConst_lev{lev+1}*RTCtoRTF{lev};
    A21AsymConst_lev{lev}=A12AsymConst_lev{lev}';
    A22AsymConst_lev{lev}=RTCtoRTF{lev}'*A22AsymConst_lev{lev+1}*RTCtoRTF{lev};
end
end


% create proper level dipendent matrices
for lev=1:L
    A11_lev{lev}=A11Eq_lev{lev}+A11AsymConst_lev{lev};
    A12_lev{lev}=A12Eq_lev{lev}+A12AsymConst_lev{lev};
    A21_lev{lev}=A21Eq_lev{lev}+A21AsymConst_lev{lev};
    A22_lev{lev}=A22Eq_lev{lev}+A22AsymConst_lev{lev};
end



% multiply rhs by h^2
b=b.*h(L)^2;


% add bc to finer level
E_remove=mesh{L}.E_remove;
E_label=mesh{L}.E_label;
A11_lev{L}(E_remove,:)=0; A12_lev{L}(E_remove,:)=0; 
A21_lev{L}(E_remove,:)=0; A22_lev{L}(E_remove,:)=0;
EcontBC=0;
NE=mesh{L}.NE;
type_of_dof=2;
for ii=E_remove
    U_xy=[0;0];
    EcontBC=EcontBC+1;
    tmp=add_boundary_bc_elasticity2D(U_xy,type_of_dof,E_label(EcontBC), contact);  
    coeff = RT_dirichlet_coeff(E_remove(EcontBC), mesh{L});
    b(ii)=tmp(1)/coeff;
    b(ii+NE)=tmp(2)/coeff;
    A11_lev{lev}(ii,ii)=1;  
    A22_lev{lev}(ii,ii)=1; 
end

Ass{L}=[A11_lev{L} A12_lev{L};
        A21_lev{L} A22_lev{L}];



% add bc to rougher levels
for lev=1:L-1
    
    NE=mesh{lev}.NE;
    E_remove=mesh{lev}.E_remove;
    
    A11_lev{lev}(E_remove,:)=0; 
    A12_lev{lev}(E_remove,:)=0; 
    A21_lev{lev}(E_remove,:)=0; 
    A22_lev{lev}(E_remove,:)=0; 
       
    for ii=E_remove
        A11_lev{lev}(ii,ii)=1;
        A22_lev{lev}(ii,ii)=1;
    end
    
    
Ass{lev}=[A11_lev{lev} A12_lev{lev};
          A21_lev{lev} A22_lev{lev}];
end
















end