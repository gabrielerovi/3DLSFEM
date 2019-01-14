function [A_lev,b_lev,C_lev,RTCtoRTF,P1CtoP1F,P1toRT]=RT_create_levels2D(A,b,C,mesh,RTCtoRTF,P1CtoP1F,P1toRT)

L=size(mesh);
L=L(1);

A11_lev=cell(L,1); A11_lev{L}=A11;
A12_lev=cell(L,1); A12_lev{L}=A12;
A22_lev=cell(L,1); A22_lev{L}=A22;
M_lev=cell(L,1);
M_lev{L}=[A11 A12;
          A12'A22];

for lev=1:L-1
  A11_lev{lev}  =RTCtoRTF{lev}' * A11_lev{lev+1} * RTCtoRTF{lev};
  A12_lev{lev}  =RTCtoRTF{lev}' * A12_lev{lev+1} * P1CtoP1F{lev};
  A22_lev{lev}  =P1CtoP1F{lev}' * A22_lev{lev+1} * P1CtoP1F{lev};
  
  M_lev{lev}=[A11_lev{lev}  A12_lev{lev};
              A12_lev{lev}' A22_lev{lev}];
          
          
    
  Z12=zeros(mesh{lev+1}.NE, mesh{lev}.N);  
  Z21=zeros(mesh{lev+1}.N, mesh{lev}.NE);  
  CtoF_lev{lev}=[RTCtoRTF{lev}    Z_P1           ;
                 Z_P1             P1CtoP1F{lev} ]; 
end




A_lev=cell(L,1);
C_lev=cell(L,1);
CtoF_lev=cell(L-1,1);
RTP1toP1P1_lev=cell(L,1);

b_dirichlet=zeros(length(A(:,1)),1);
A_lev{L}=A;
% create projection matrix C_lev =[ C 0;
%                                   0 I;] : RT + P1 -> ND(P1) + P1

Z_P1=zeros(mesh{L}.N,mesh{L}.N);
I_P1=eye(mesh{L}.N,mesh{L}.N);

% create main matrix A_lev: (RT+P1) -> (RT+P1)^*
A_lev{L}=[A11  A12;
          A12' A22];
% create bilinear matrix C_lev: (ND/P1 + P1) -> (ND/P1 + P1)^*, C =[ C 0;
%                                                                    0 I;] 
C_lev{L}=[C     Z_P1;
          Z_P1  I_P1   ];

% create coarse to fine projection matrix CtoF_lev: (RT_C + P1_C) -> (RT_F + P1_F), C =[ P1toRT 0;
%                                                                                        0      I;] 
for lev=1:L-1
  Z12=zeros(mesh{lev+1}.NE, mesh{lev}.N);  
  Z21=zeros(mesh{lev+1}.N, mesh{lev}.NE);  
  CtoF_lev{lev}=[RTCtoRTF{lev}    Z_P1           ;
                 Z_P1             P1CtoP1F{lev} ]; 
end
% create projection matrix RTP1toP1P1_lev: (ND/P1 + P1) -> (RT + P1), C =[ P1toRT 0;
%                                                                          0      I;] 
for lev=1:L
  Z=zeros(mesh{lev}.NE, mesh{lev}.N);  
  I=eye(mesh{lev}.N, mesh{lev}.N);  
  RTP1toP1P1_lev{lev}=[P1toRT{lev}    Z ;
                       Z'             I ]; 
end     
  
  
  
if(L>1)
for lev=L-1:-1:1   
    A_lev{lev}=RTCtoRTF{lev}'*A_lev{lev+1}*RTCtoRTF{lev};          
    C_lev{lev}=[P1CtoP1F{lev}'*C_lev{lev+1}*P1CtoP1F{lev}];
end
end


    Z_P1=zeros(mesh{lev}.N,mesh{lev}.N);
    C_lev{L}=[C_lev{lev}         Z_P1;
              Z_P1               I_P1];


E_bc=mesh{L}.E_bc;
E_dirichlet=mesh{L}.E_dirichlet;
E_remove=mesh{L}.E_remove;
N_remove=mesh{L}.N_remove;
internal=0;
boundary=0;
type_of_dof=2;

for eF=1:mesh{L}.NE
    if(E_dirichlet(eF)>0)
        b_dirichlet(eF)=dirichlet_bc(E_bc(eF),type_of_dof); 
    else
% for each internal dof, compute A_{i,:)b_dirichlet(i), that is known and
% can be put in the rhs
        internal=internal+1;
    end
end
rhs=zeros(internal,1);
internal=0;
type_of_dof=2;
for eF=1:mesh{L}.NE
    if(E_dirichlet(eF)==0 )
% for each internal dof, compute A_{i,:)b_dirichlet(i), that is known and
% can be put in the rhs
        internal=internal+1;
        rhs(internal)=A(eF,:)*b_dirichlet;
    end
end


A_lev{L}(E_remove,:)=[];
A_lev{L}(:,E_remove)=[];
b(E_remove)=[];
b_lev=b-rhs;

C_lev{L}(N_remove,:)=[];
C_lev{L}(:,N_remove)=[];

P1toRT{L}(E_remove,:)=[];
P1toRT{L}(:,N_remove)=[];

for lev=1:L-1
    NC_remove=mesh{lev}.N_remove; 
    EC_remove=mesh{lev}.E_remove;
    NF_remove=mesh{lev+1}.N_remove;
    EF_remove=mesh{lev+1}.E_remove;
    
    C_lev{lev}(NC_remove,:)=[];
    C_lev{lev}(:,NC_remove)=[];
    A_lev{lev}(EC_remove,:)=[];
    A_lev{lev}(:,EC_remove)=[];
    
    RTCtoRTF{lev}(EF_remove,:)=[];
    RTCtoRTF{lev}(:,EC_remove)=[];
    
    P1CtoP1F{lev}(NF_remove,:)=[];
    P1CtoP1F{lev}(:,NC_remove)=[];
    
    P1toRT{lev}(EC_remove,:)=[];
    P1toRT{lev}(:,NC_remove)=[];
    
end



end