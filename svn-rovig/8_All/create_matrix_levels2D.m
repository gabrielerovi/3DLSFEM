function [M_tot,b_tot,P1CtoP1F_block,RTCtoRTF_block,NDtoRT_block,A11_block,A12_block,A22_block,b1_block,b2_block,C_block,...
    A11_lev,A12_lev,A13_lev,A14_lev,A22_lev,A23_lev,A24_lev,A33_lev,A34_lev,A44_lev,...
          b1_lev,b2_lev,b3_lev,b4_lev,C11_lev,C12_lev,C22_lev]=create_matrix_levels2D ...
          (A11,A12,A13,A14,A22,A23,A24,A33,A34,A44,C11,C12,C22,b1,b2,b3,b4,mesh,RTCtoRTF,P1CtoP1F,NDCtoNDF,NDtoRT)

L=size(mesh);
L=L(1);
A11_lev=cell(L,1);
A12_lev=cell(L,1);
A13_lev=cell(L,1);
A14_lev=cell(L,1);

A22_lev=cell(L,1);
A23_lev=cell(L,1);
A24_lev=cell(L,1);

A33_lev=cell(L,1);
A34_lev=cell(L,1);
A44_lev=cell(L,1);

A11_block=cell(L,1);
A12_block=cell(L,1);
A22_block=cell(L,1);

C_block=cell(L,1);

P1CtoP1F_block=cell(L-1,1);
RTCtoRTF_block=cell(L-1,1);            
NDtoRT_block=cell(L,1);
    
    
C11_lev=cell(L,1);
C12_lev=cell(L,1);
C22_lev=cell(L,1);

A11_lev{L}=A11;
A12_lev{L}=A12;
A13_lev{L}=A13;
A14_lev{L}=A14;

A22_lev{L}=A22;
A23_lev{L}=A23;
A24_lev{L}=A24;

A33_lev{L}=A33;
A34_lev{L}=A34;
A44_lev{L}=A44;

C11_lev{L}=C11;
C12_lev{L}=C12;
C22_lev{L}=C22;

if(L>1)
for lev=L-1:-1:1   
    A11_lev{lev}=RTCtoRTF{lev}'*A11_lev{lev+1}*RTCtoRTF{lev};
    A12_lev{lev}=RTCtoRTF{lev}'*A12_lev{lev+1}*RTCtoRTF{lev};
    A13_lev{lev}=RTCtoRTF{lev}'*A13_lev{lev+1}*P1CtoP1F{lev};   
    A14_lev{lev}=RTCtoRTF{lev}'*A14_lev{lev+1}*P1CtoP1F{lev};
    
    A22_lev{lev}=RTCtoRTF{lev}'*A22_lev{lev+1}*RTCtoRTF{lev};
    A23_lev{lev}=RTCtoRTF{lev}'*A23_lev{lev+1}*P1CtoP1F{lev};   
    A24_lev{lev}=RTCtoRTF{lev}'*A24_lev{lev+1}*P1CtoP1F{lev};    

    A33_lev{lev}=P1CtoP1F{lev}'*A33_lev{lev+1}*P1CtoP1F{lev}; 
    A34_lev{lev}=P1CtoP1F{lev}'*A34_lev{lev+1}*P1CtoP1F{lev}; 
    A44_lev{lev}=P1CtoP1F{lev}'*A44_lev{lev+1}*P1CtoP1F{lev}; 
    
    C11_lev{lev}=NDCtoNDF{lev}'*C11_lev{lev+1}*NDCtoNDF{lev};
    C12_lev{lev}=NDCtoNDF{lev}'*C12_lev{lev+1}*NDCtoNDF{lev};
    C22_lev{lev}=NDCtoNDF{lev}'*C22_lev{lev+1}*NDCtoNDF{lev};
    
end
end









% build the rhs with only dirichlet (edge and nodal)
% if the dirichlet (not necessarily the boundary) are on the objects edges
% 1,2, where the mesh has edges 1,2,3,4,5, then:
% E_dirichlet=[alpha1 alpha2 0 0 0]
type_of_dof=1;
[N_dirichlet,N_bc_dirichlet]=dirichlet_bc_vector(mesh{L}.N_bc,type_of_dof,mesh{L});
type_of_dof=2;
[E_dirichlet,E_bc_dirichlet]=dirichlet_bc_vector(mesh{L}.E_bc,type_of_dof,mesh{L});

% build the whole rhs, subtracting the dirichlet condition
b1_lev= b1 -  A11 * E_dirichlet(:,1)  - A12 * E_dirichlet(:,2) - A13 * N_dirichlet(:,1)  - A14 * N_dirichlet(:,2);
b2_lev= b2 -  A12'* E_dirichlet(:,1)  - A22 * E_dirichlet(:,2) - A23 * N_dirichlet(:,1)  - A24 * N_dirichlet(:,2);
b3_lev= b3 -  A13'* E_dirichlet(:,1)  - A23'* E_dirichlet(:,2) - A33 * N_dirichlet(:,1)  - A34 * N_dirichlet(:,2);
b4_lev= b4 -  A14'* E_dirichlet(:,1)  - A24'* E_dirichlet(:,2) - A34'* N_dirichlet(:,1)  - A44 * N_dirichlet(:,2);

% consider the not-dirichlet components of the rhs and remove the others
E_remove=mesh{L}.E_remove;
N_remove=mesh{L}.N_remove;

b1_lev(E_remove)=[];
b2_lev(E_remove)=[];
b3_lev(N_remove)=[];
b4_lev(N_remove)=[];


E_removeC=[];
N_removeC=[];
% for each level 
for lev=1:L
% consider the not-dirichlet components of the the matrices 
type_of_dof=2;
% we consider the nodes, but a bc==dirichlet only if it is for RT
% for this reason we use N_bc, but type_of_dof=2
[Ned_bool,Ned_Marker]=is_surface_dirichlet(mesh{lev}.N_bc,type_of_dof);
Ned_remove=find(Ned_bool);
E_remove=mesh{lev}.E_remove;
N_remove=mesh{lev}.N_remove;

% E_rem, E_remcA11_lev{lev}(E_remove,:)=[];
A11_lev{lev}(:,E_remove)=[];
A11_lev{lev}(E_remove,:)=[];

A12_lev{lev}(E_remove,:)=[];
A12_lev{lev}(:,E_remove)=[];

A22_lev{lev}(E_remove,:)=[];
A22_lev{lev}(:,E_remove)=[];

% E_rem, N_rem
A13_lev{lev}(E_remove,:)=[];
A13_lev{lev}(:,N_remove)=[];

A14_lev{lev}(E_remove,:)=[];
A14_lev{lev}(:,N_remove)=[];

A23_lev{lev}(E_remove,:)=[];
A23_lev{lev}(:,N_remove)=[];

A24_lev{lev}(E_remove,:)=[];
A24_lev{lev}(:,N_remove)=[];

% N_rem, N_rem
A33_lev{lev}(N_remove,:)=[];
A33_lev{lev}(:,N_remove)=[];
 
A34_lev{lev}(N_remove,:)=[];
A34_lev{lev}(:,N_remove)=[];

A44_lev{lev}(N_remove,:)=[];
A44_lev{lev}(:,N_remove)=[];

A11_block{lev}=[A11_lev{lev}  A12_lev{lev};
                A12_lev{lev}' A22_lev{lev}];
A12_block{lev}=[A13_lev{lev}  A14_lev{lev};
                A23_lev{lev} A24_lev{lev}];
A22_block{lev}=[A33_lev{lev}  A34_lev{lev};
                A34_lev{lev}' A44_lev{lev}];
C_block{lev}=[C11_lev{lev}  C12_lev{lev};
              C12_lev{lev}' C22_lev{lev};];

NDtoRT{lev}(E_remove,:)=[];
ss=size(NDtoRT{lev});                                      
NDtoRT_block{lev}=[NDtoRT{lev}        zeros(ss(1),ss(2));
                   zeros(ss(1),ss(2)) NDtoRT{lev}  ];
                   
if(lev>1)          
P1CtoP1F{lev-1}(N_remove,:)=[];     
P1CtoP1F{lev-1}(:,N_removeC)=[];
RTCtoRTF{lev-1}(E_remove,:)=[];
RTCtoRTF{lev-1}(:,E_removeC)=[];
                 
    ss=size(P1CtoP1F{lev-1});
    
    P1CtoP1F_block{lev-1}=[P1CtoP1F{lev-1},      zeros(ss(1),ss(2));
                         zeros(ss(1),ss(2)), P1CtoP1F{lev-1};];
    ss=size(RTCtoRTF{lev-1});
    RTCtoRTF_block{lev-1}=[RTCtoRTF{lev-1},      zeros(ss(1),ss(2));
                         zeros(ss(1),ss(2)), RTCtoRTF{lev-1}];
end
      
E_removeC=E_remove;
N_removeC=N_remove;
end




b1_block=[b1_lev;b2_lev];
b2_block=[b3_lev;b4_lev];



M_tot=[A11, A12, A13, A14;
       A12',A22, A23, A24;
       A13',A23',A33, A34;
       A14',A24',A34',A44;];
 
b_tot=[b1;b2;b3;b4];
NE=mesh{L}.NE;
N=mesh{L}.N;
type_of_dof=2;
for ee=1:NE
    [bool_loc,Marker_loc]=is_surface_dirichlet(mesh{L}.E_bc(ee),type_of_dof);
    if(bool_loc>0)
        ee2=ee+NE;
        M_tot(ee,:)=0;     
        M_tot(ee,ee)=1; 
        b_tot(ee)=E_dirichlet(ee,1); 
        M_tot(ee2,:)=0;  
        M_tot(ee2,ee2)=1;
        b_tot(ee2)=E_dirichlet(ee,2);
        
    end
end

type_of_dof=1;
for nn=1:N
    [bool_loc,Marker_loc]=is_surface_dirichlet(mesh{L}.N_bc(nn),type_of_dof);
    if(bool_loc>0)
        nn1=nn+2*NE;
        nn2=nn1+N;
        M_tot(nn1,:)=0;     
        M_tot(nn1,nn1)=1; 
        b_tot(nn1)=N_dirichlet(nn,1); 
        M_tot(nn2,:)=0;  
        M_tot(nn2,nn2)=1;
        b_tot(nn2)=N_dirichlet(nn,2);
        
    end
end

end