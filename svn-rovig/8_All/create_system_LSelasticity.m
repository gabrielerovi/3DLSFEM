
function [A,b,AFine,P, Ant,bnt,Antbc,AFinenobc,AFinenobcnt,bnt_lev,C]=create_system_LSelasticity(parameters,mesh,h,P)

% problem coefficients
input_name=parameters.input_name;
qrule=parameters.qrule;
L=length(mesh);
NE=mesh{L}.NE;
N=mesh{L}.N;
N_components=parameters.N_components;
E_components=parameters.E_components;
E_remove=mesh{L}.E_remove;
E_label=mesh{L}.E_label;
N_remove=mesh{L}.N_remove;
N_label=mesh{L}.N_label;
dim_problem=N_components+E_components;

% coefficients 
C_eq=parameters.C_eq;
C_const=parameters.C_const;
C_asym=parameters.C_asym;
qrule=parameters.qrule;
alpha=parameters.alpha;
beta=parameters.beta;
contact=parameters.contact;
penalty=parameters.penalty;

%projection
RTCtoRTF= P.RTCtoRTF;
P1CtoP1F= P.P1CtoP1F;


if(parameters.h_dependence==false)
    h=ones(L,1);
end

if(parameters.lump_asym==true &&parameters.lump_diag_asym ==false)
    C_lump=C_asym;
    C_lump_diag=0;
    C_asym_use_diag=0;
    C_asym_use_extradiag=0;
elseif(parameters.lump_asym==false &&parameters.lump_diag_asym ==true)
    C_lump_diag=C_asym;
    C_lump=0;
    C_asym_use_diag=0; 
    C_asym_use_extradiag=0;
else
    C_asym_use_diag=1;
    C_lump=0;
    C_lump_diag=0;
    C_asym_use_extradiag=1;
end


switch N_components
    case 0
        N_remove=[];
        N_label=[];
    otherwise
        N_remove=[mesh{L}.N_remove];
        N_label=mesh{L}.N_label;
end

switch E_components
    case 0
        E_remove=[];
        E_label=[];
    otherwise
        E_remove=[mesh{L}.E_remove];
        E_label=mesh{L}.E_label;
end

Aeq=cell(L,dim_problem,dim_problem);
Aasym=cell(L,dim_problem,dim_problem);
A=cell(L,dim_problem,dim_problem);

% stress block part
% we consider it separately because the equilibrium term can be h-deepndent 
Aeq{L,1,1}=h(L)^2 * assembling2DRTRT(mesh{L},qrule,alpha,beta,1,1,C_eq,0,0,input_name);
Aeq{L,1,2}=h(L)^2 *assembling2DRTRT(mesh{L},qrule,alpha,beta,1,2,C_eq,0,0,input_name);
Aeq{L,2,1}=Aeq{L,1,2}';
Aeq{L,2,2}=h(L)^2 *assembling2DRTRT(mesh{L},qrule,alpha,beta,2,2,C_eq,0,0,input_name);

Aasym{L,1,1}=assembling2DRTRT(mesh{L},qrule,alpha,beta,1,1,0,C_const,C_asym,input_name);
Aasym{L,1,2}=assembling2DRTRT(mesh{L},qrule,alpha,beta,1,2,0,C_const,C_asym,input_name);
Aasym{L,2,1}=Aasym{L,1,2}';
Aasym{L,2,2}=assembling2DRTRT(mesh{L},qrule,alpha,beta,2,2,0,C_const,C_asym,input_name);

[DiagAsym1,DiagAsym2]=Assembling_Global_Lumped_Asymmetry_Sigma(mesh{L},qrule,C_asym);

Aasym{L,1,1}= C_asym_use_diag * Aasym{L,1,1} + C_lump_diag * diag(diag(Aasym{L,1,1})) + C_lump * DiagAsym1;
Aasym{L,1,2}= C_asym_use_extradiag * Aasym{L,1,2} ;
Aasym{L,2,1}= C_asym_use_extradiag * Aasym{L,2,1} ;
Aasym{L,2,2}= C_asym_use_diag * Aasym{L,2,2} + C_lump_diag * diag(diag(Aasym{L,2,2})) + C_lump * DiagAsym2;




A{L,1,3}=assembling2DRTP1(mesh{L},qrule,alpha,beta,1,3,C_eq,C_const,input_name);
A{L,1,4}=assembling2DRTP1(mesh{L},qrule,alpha,beta,1,4,C_eq,C_const,input_name);
A{L,2,3}=assembling2DRTP1(mesh{L},qrule,alpha,beta,2,3,C_eq,C_const,input_name);
A{L,2,4}=assembling2DRTP1(mesh{L},qrule,alpha,beta,2,4,C_eq,C_const,input_name);

A{L,3,1}=assembling2DRTP1(mesh{L},qrule,alpha,beta,3,1,C_eq,C_const,input_name);
A{L,3,2}=assembling2DRTP1(mesh{L},qrule,alpha,beta,3,2,C_eq,C_const,input_name);
A{L,3,3}=assembling2DP1P1(mesh{L},qrule,alpha,beta,3,3,C_eq,C_const,input_name);
A{L,3,4}=assembling2DP1P1(mesh{L},qrule,alpha,beta,3,4,C_eq,C_const,input_name);

A{L,4,1}=assembling2DRTP1(mesh{L},qrule,alpha,beta,4,1,C_eq,C_const,input_name);
A{L,4,2}=assembling2DRTP1(mesh{L},qrule,alpha,beta,4,2,C_eq,C_const,input_name);
A{L,4,3}=A{L,3,4}';
A{L,4,4}=assembling2DP1P1(mesh{L},qrule,alpha,beta,4,4,C_eq,C_const,input_name);

b1=h(L)^2 * assemblingb11(mesh{L},qrule,parameters.force1,C_eq,0);
b2=h(L)^2 * assemblingb11(mesh{L},qrule,parameters.force2,C_eq,0);
b3=assemblingb22(mesh{L},qrule,parameters.fzero,C_eq,C_const);
b4=assemblingb22(mesh{L},qrule,parameters.fzero,C_eq,C_const);

% if(contact==1)
%     
%     [A,b1,b2] =...
%        ComplementarityCondition(A,L,b1,b2,mesh{L},parameters);
%     
% end

b=[b1;b2;b3;b4];




if(L>1)




for lev=L-1:-1:1   
    
    % stress block 11,12,21,22
    Aeq{lev,1,1}=(h(lev)/h(lev+1))^2.* RTCtoRTF{lev}'*Aeq{lev+1,1,1}*RTCtoRTF{lev};
    Aeq{lev,1,2}=(h(lev)/h(lev+1))^2.* RTCtoRTF{lev}'*Aeq{lev+1,1,2}*RTCtoRTF{lev};
    Aeq{lev,2,1}=Aeq{lev,1,2}';
    Aeq{lev,2,2}=(h(lev)/h(lev+1))^2.* RTCtoRTF{lev}'*Aeq{lev+1,2,2}*RTCtoRTF{lev};
 
    Aasym{lev,1,1}= RTCtoRTF{lev}'*Aasym{lev+1,1,1}*RTCtoRTF{lev};
    Aasym{lev,1,2}= RTCtoRTF{lev}'*Aasym{lev+1,1,2}*RTCtoRTF{lev};
    Aasym{lev,2,1}=Aasym{lev,1,2}';
    Aasym{lev,2,2}= RTCtoRTF{lev}'*Aasym{lev+1,2,2}*RTCtoRTF{lev};
    
    
    % continue with the rest of the matrix...
    A{lev,1,3}=RTCtoRTF{lev}'*A{lev+1,1,3}*P1CtoP1F{lev};   
    A{lev,1,4}=RTCtoRTF{lev}'*A{lev+1,1,4}*P1CtoP1F{lev};      
    A{lev,2,3}=RTCtoRTF{lev}'*A{lev+1,2,3}*P1CtoP1F{lev};   
    A{lev,2,4}=RTCtoRTF{lev}'*A{lev+1,2,4}*P1CtoP1F{lev};  
    
    
    A{lev,3,1}=P1CtoP1F{lev}'*A{lev+1,3,1}*RTCtoRTF{lev};
    A{lev,3,2}=P1CtoP1F{lev}'*A{lev+1,3,2}*RTCtoRTF{lev};
    A{lev,3,3}=P1CtoP1F{lev}'*A{lev+1,3,3}*P1CtoP1F{lev}; 
    A{lev,3,4}=P1CtoP1F{lev}'*A{lev+1,3,4}*P1CtoP1F{lev}; 
    
    
    A{lev,4,1}=P1CtoP1F{lev}'*A{lev+1,4,1}*RTCtoRTF{lev};
    A{lev,4,2}=P1CtoP1F{lev}'*A{lev+1,4,2}*RTCtoRTF{lev};
    A{lev,4,3}=A{lev,3,4}';
    A{lev,4,4}=P1CtoP1F{lev}'*A{lev+1,4,4}*P1CtoP1F{lev}; 

    
end
end




for lev=L
    A{lev,1,1}=Aeq{lev,1,1} + Aasym{lev,1,1};
    A{lev,1,2}=Aeq{lev,1,2} + Aasym{lev,1,2};
    A{lev,2,1}=Aeq{lev,2,1} + Aasym{lev,2,1};
    A{lev,2,2}=Aeq{lev,2,2} + Aasym{lev,2,2};
% in case of contact, add proper terms from the complementarity condition
% < sigma n n, u n - g> = 0 added to the functional


if(penalty==1)
    
    [A,b3,b4] =...
       ContactPenalty(A,L,b3,b4,mesh{L},parameters);
    
end

end

for lev=1:L-1
    A{lev,1,1}=Aeq{lev,1,1} + Aasym{lev,1,1};
    A{lev,1,2}=Aeq{lev,1,2} + Aasym{lev,1,2};
    A{lev,2,1}=Aeq{lev,2,1} + Aasym{lev,2,1};
    A{lev,2,2}=Aeq{lev,2,2} + Aasym{lev,2,2};
% in case of contact, add proper terms from the complementarity condition
% < sigma n n, u n - g> = 0 added to the functional
% if(contact==1)
%     
%     [A,b1,b2] =...
%        ComplementarityCondition(A,L,b1,b2,mesh{lev},parameters);
%     
% end
end

lev=L;
AFinenobc=[A{lev,1,1} A{lev,1,2} A{lev,1,3} A{lev,1,4};
              A{lev,2,1} A{lev,2,2} A{lev,2,3} A{lev,2,4};
              A{lev,3,1} A{lev,3,2} A{lev,3,3} A{lev,3,4};
              A{lev,4,1} A{lev,4,2} A{lev,4,3} A{lev,4,4};];
          
          
householder=1;
for lev=1:L
Antnobc{lev}=[A{lev,1,1} A{lev,1,2} A{lev,1,3} A{lev,1,4};
              A{lev,2,1} A{lev,2,2} A{lev,2,3} A{lev,2,4};
              A{lev,3,1} A{lev,3,2} A{lev,3,3} A{lev,3,4};
              A{lev,4,1} A{lev,4,2} A{lev,4,3} A{lev,4,4};];
          
     [M_Normal_Tangent,M_Normal_TangentT] = MatrixOnGammaCwithNormalTangentComponents(mesh{lev},householder);
     Antnobc{lev} = M_Normal_Tangent * Antnobc{lev} * M_Normal_Tangent';
     LLL=[0, mesh{lev}.NE, 2*mesh{lev}.NE, mesh{lev}.N + 2 * mesh{lev}.NE, 2 * mesh{lev}.N + 2 * mesh{lev}.NE ];
     
     for mm=1:4
     for nn=1:4
         Ant{lev,mm,nn}=Antnobc{lev}(LLL(mm)+1:LLL(mm+1),LLL(nn)+1:LLL(nn+1));
         Antbc{lev,mm,nn}=Ant{lev,mm,nn};
     end
     end
     
     if(contact==1)
    
        [C{lev}.Ant_complementarity,C{lev}.bnt_complementarity] =...
           ComplementarityConditionNT(lev,mesh,parameters);
     end
end

AFinenobcnt= [Ant{L,1,1}, Ant{L,1,2}, Ant{L,1,3}, Ant{L,1,4};
              Ant{L,2,1}, Ant{L,2,2}, Ant{L,2,3}, Ant{L,2,4};
              Ant{L,3,1}, Ant{L,3,2}, Ant{L,3,3}, Ant{L,3,4};
              Ant{L,4,1}, Ant{L,4,2}, Ant{L,4,3}, Ant{L,4,4};];
        
     
     
     
     
% add bc to finer level
lev=L;
A{L,1,1}(E_remove,:)=0; 
A{L,1,2}(E_remove,:)=0; 
A{L,1,3}(E_remove,:)=0; 
A{L,1,4}(E_remove,:)=0; 

A{L,2,1}(E_remove,:)=0; 
A{L,2,2}(E_remove,:)=0;
A{L,2,3}(E_remove,:)=0;
A{L,2,4}(E_remove,:)=0;

A{L,3,1}(N_remove,:)=0; 
A{L,3,2}(N_remove,:)=0;
A{L,3,3}(N_remove,:)=0;
A{L,3,4}(N_remove,:)=0;

A{L,4,1}(N_remove,:)=0; 
A{L,4,2}(N_remove,:)=0;
A{L,4,3}(N_remove,:)=0;
A{L,4,4}(N_remove,:)=0;

EcontBC=0;
type_of_dof=2;
for ii=E_remove
    U_xy=[0;0];
    EcontBC=EcontBC+1;
    tmp=add_boundary_bc_elasticity2D(U_xy,type_of_dof,E_label(EcontBC), 0);  
    coeff = RT_dirichlet_coeff(E_remove(EcontBC), mesh{L});
    jj1=ii;
    jj2=ii+NE;    
    b(jj1)=tmp(1)/coeff;
    b(jj2)=tmp(2)/coeff;
    A{lev,1,1}(ii,ii)=1;  
    A{lev,2,2}(ii,ii)=1; 
end

NcontBC=0;
type_of_dof=1;
for ii=N_remove
    U_xy=[0;0];
    NcontBC=NcontBC+1;
    tmp=add_boundary_bc_elasticity2D(U_xy,type_of_dof,N_label(NcontBC), 0);   
    jj1=ii+2*NE;
    jj2=ii+2*NE+N;
    b(jj1)=tmp(1);
    b(jj2)=tmp(2);
    A{lev,3,3}(ii,ii)=1;  
    A{lev,4,4}(ii,ii)=1; 

end


% add bc to rougher levels
for lev=1:L-1

E_remove=mesh{lev}.E_remove;
N_remove=mesh{lev}.N_remove;
A{lev,1,1}(E_remove,:)=0; 
A{lev,1,2}(E_remove,:)=0; 
A{lev,1,3}(E_remove,:)=0; 
A{lev,1,4}(E_remove,:)=0; 

A{lev,2,1}(E_remove,:)=0; 
A{lev,2,2}(E_remove,:)=0;
A{lev,2,3}(E_remove,:)=0;
A{lev,2,4}(E_remove,:)=0;

A{lev,3,1}(N_remove,:)=0; 
A{lev,3,2}(N_remove,:)=0;
A{lev,3,3}(N_remove,:)=0;
A{lev,3,4}(N_remove,:)=0;

A{lev,4,1}(N_remove,:)=0; 
A{lev,4,2}(N_remove,:)=0;
A{lev,4,3}(N_remove,:)=0;
A{lev,4,4}(N_remove,:)=0; 
       
    for ii=E_remove 
        A{lev,1,1}(ii,ii)=1;
        A{lev,2,2}(ii,ii)=1;
    end
    for ii=N_remove
        A{lev,3,3}(ii,ii)=1;  
        A{lev,4,4}(ii,ii)=1; 
    end
    
    
end



  
AFine=[A{L,1,1} A{L,1,2} A{L,1,3} A{L,1,4};
       A{L,2,1} A{L,2,2} A{L,2,3} A{L,2,4};
       A{L,3,1} A{L,3,2} A{L,3,3} A{L,3,4};
       A{L,4,1} A{L,4,2} A{L,4,3} A{L,4,4};];




for lev=1:L
    E_remove=mesh{lev}.E_remove;
    N_remove=mesh{lev}.N_remove;
    E_contact=mesh{lev}.E_contact;
    
         for ii=1:2
             for jj=1:4
                 Antbc{lev,ii,jj}(E_remove,:)=0;
                 Antbc{lev,ii+2,jj}(N_remove,:)=0;
             end
         end
    
         for ee1=E_remove
             for ii=1:2
             Antbc{lev,ii,ii}(ee1,ee1)=1;
             end
         end
       
         
         for nn1=N_remove
             for ii=3:4
             Antbc{lev,ii,ii}(nn1,nn1)=1;
             end
         end  
         
         for ee1=E_contact
             for mm=1:4
%               % we remove also the normal component, like we were computing pure dirichlet problem  
%              Antbc{lev,1,mm}(ee1,:)=0;
%            
% since friction ==0, we can just make the matrux symmetric
             Antbc{lev,2,mm}(ee1,:)=0;
             Antbc{lev,mm,2}(:,ee1)=0;
             end
%              
%              Antbc{lev,1,1}(ee1,ee1)=1;
             Antbc{lev,2,2}(ee1,ee1)=1;
         end
         
         
  
         
         

end






for lev=1:L
% here we compute the rhs b for each mesh
b_lev1=h(lev)^2 * assemblingb11(mesh{lev},qrule,parameters.force1,C_eq,0);
b_lev2=h(lev)^2 * assemblingb11(mesh{lev},qrule,parameters.force2,C_eq,0);
b_lev3=assemblingb22(mesh{lev},qrule,parameters.fzero,C_eq,C_const);
b_lev4=assemblingb22(mesh{lev},qrule,parameters.fzero,C_eq,C_const);
b_lev{lev}=[b_lev1;b_lev2;b_lev3;b_lev4];

E_remove=mesh{lev}.E_remove;
N_remove=mesh{lev}.N_remove;
E_label=mesh{lev}.E_label;
N_label=mesh{lev}.N_label;
NE=mesh{lev}.NE;
N=mesh{lev}.N;

EcontBC=0;
type_of_dof=2;
for ii=E_remove
    U_xy=[0;0];
    EcontBC=EcontBC+1;
    tmp=add_boundary_bc_elasticity2D(U_xy,type_of_dof,E_label(EcontBC), 0);  
    coeff = RT_dirichlet_coeff(E_remove(EcontBC), mesh{lev});
    jj1=ii;
    jj2=ii+NE;    
    b_lev{lev}(jj1)=tmp(1)/coeff;
    b_lev{lev}(jj2)=tmp(2)/coeff;
end

NcontBC=0;
type_of_dof=1;
for ii=N_remove
    U_xy=[0;0];
    NcontBC=NcontBC+1;
    tmp=add_boundary_bc_elasticity2D(U_xy,type_of_dof,N_label(NcontBC), 0);   
    jj1=ii+2*NE;
    jj2=ii+2*NE+N;
    b_lev{lev}(jj1)=tmp(1);
    b_lev{lev}(jj2)=tmp(2);
end

[M_Normal_Tangent,M_Normal_TangentT] = MatrixOnGammaCwithNormalTangentComponents(mesh{lev},householder);
bnt_lev{lev}=M_Normal_Tangent*b_lev{lev};

end
bnt=bnt_lev{L};
end



 
    






