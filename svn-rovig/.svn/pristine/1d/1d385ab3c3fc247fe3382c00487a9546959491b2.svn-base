
function [A,b,AFine,P]=create_system_LSstressblocks(parameters,mesh,h,P)

% problem coefficients
input_name=parameters.input_name;
qrule=parameters.qrule;
L=length(mesh);
N_components=parameters.N_components;
E_components=parameters.E_components;
E_remove=mesh{L}.E_remove;
E_label=mesh{L}.E_label;
dim_problem=N_components+E_components;

% coefficients 
C_eq=parameters.C_eq;
C_const=parameters.C_const;
C_asym=parameters.C_asym;
qrule=parameters.qrule;
alpha=parameters.alpha;
beta=parameters.beta;
contact=parameters.contact;

%projection
RTCtoRTF=P.RTCtoRTF;



if(parameters.h_dependence==false)
    h=ones(L,1);
end

if(parameters.lump_asym==true &&parameters.lump_diag_asym ==false)
    C_lump=C_asym;
    C_lump_diag=0;
    C_asym=0;
elseif(parameters.lump_asym==true &&parameters.lump_diag_asym ==true)
    C_lump_diag=C_asym;
    C_lump=0;
    C_asym=0;    
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

Aeq{L,1,1}=h(L)^2 * assembling2DRTRT(mesh{L},qrule,alpha,beta,1,1,C_eq,0,0,input_name);
Aeq{L,1,2}=h(L)^2 *assembling2DRTRT(mesh{L},qrule,alpha,beta,1,2,C_eq,0,0,input_name);
Aeq{L,2,1}=Aeq{L,1,2}';
Aeq{L,2,2}=h(L)^2 *assembling2DRTRT(mesh{L},qrule,alpha,beta,2,2,C_eq,0,0,input_name);

Aasym{L,1,1}=assembling2DRTRT(mesh{L},qrule,alpha,beta,1,1,0,C_const,C_asym,input_name);
Aasym{L,1,2}=assembling2DRTRT(mesh{L},qrule,alpha,beta,1,2,0,C_const,C_asym,input_name);
Aasym{L,2,1}=Aasym{L,1,2}';
Aasym{L,2,2}=assembling2DRTRT(mesh{L},qrule,alpha,beta,2,2,0,C_const,C_asym,input_name);

[DiagAsym1,DiagAsym2]=Assembling_Global_Lumped_Asymmetry_Sigma(mesh{L},qrule,C_asym);

Aasym{L,1,1}=Aasym{L,1,1} + DiagAsym1 + C_lump_diag * diag(diag(Aasym{L,1,1}));
Aasym{L,2,2}=Aasym{L,2,2} + DiagAsym2 + C_lump_diag * diag(diag(Aasym{L,2,2}));



b1=h(L)^2 * assemblingb11(mesh{L},qrule,parameters.force1,C_eq,0);
b2=h(L)^2 * assemblingb11(mesh{L},qrule,parameters.force2,C_eq,0);

b=[b1;b2];

if(L>1)
for lev=L-1:-1:1   
    
    Aeq{lev,1,1}=(h(lev)/h(lev+1))^2.* RTCtoRTF{lev}'*Aeq{lev+1,1,1}*RTCtoRTF{lev};
    Aeq{lev,1,2}=(h(lev)/h(lev+1))^2.* RTCtoRTF{lev}'*Aeq{lev+1,1,2}*RTCtoRTF{lev};
    Aeq{lev,2,1}=Aeq{lev,1,2}';
    Aeq{lev,2,2}=(h(lev)/h(lev+1))^2.* RTCtoRTF{lev}'*Aeq{lev+1,2,2}*RTCtoRTF{lev};
 
    Aasym{lev,1,1}= RTCtoRTF{lev}'*Aasym{lev+1,1,1}*RTCtoRTF{lev};
    Aasym{lev,1,2}= RTCtoRTF{lev}'*Aasym{lev+1,1,2}*RTCtoRTF{lev};
    Aasym{lev,2,1}=Aasym{lev,1,2}';
    Aasym{lev,2,2}= RTCtoRTF{lev}'*Aasym{lev+1,2,2}*RTCtoRTF{lev};
end
end

for lev=1:L
    A{lev,1,1}=Aeq{lev,1,1} + Aasym{lev,1,1};
    A{lev,1,2}=Aeq{lev,1,2} + Aasym{lev,1,2};
    A{lev,2,1}=Aeq{lev,2,1} + Aasym{lev,2,1};
    A{lev,2,2}=Aeq{lev,2,2} + Aasym{lev,2,2};
end




% add bc to finer level
A{L,1,1}(E_remove,:)=0; A{L,1,2}(E_remove,:)=0; 
A{L,2,1}(E_remove,:)=0; A{L,2,2}(E_remove,:)=0;
EcontBC=0;
NE=mesh{L}.NE;
type_of_dof=2;
for ii=E_remove
    U_xy=[0;0];
    EcontBC=EcontBC+1;
    tmp=add_boundary_bc_elasticity2D(U_xy,type_of_dof,E_label(EcontBC), 0);  
    coeff = RT_dirichlet_coeff(E_remove(EcontBC), mesh{L});
    b(ii)=tmp(1)/coeff;
    b(ii+NE)=tmp(2)/coeff;
    A{lev,1,1}(ii,ii)=1;  
    A{lev,2,2}(ii,ii)=1; 
end

% add bc to rougher levels
for lev=1:L-1
    
    NE=mesh{lev}.NE;
    E_remove=mesh{lev}.E_remove;
    
    A{lev,1,1}(E_remove,:)=0; 
    A{lev,1,2}(E_remove,:)=0; 
    A{lev,2,1}(E_remove,:)=0; 
    A{lev,2,2}(E_remove,:)=0; 
       
    for ii=E_remove
        A{lev,1,1}(ii,ii)=1;
        A{lev,2,2}(ii,ii)=1;
    end
    
end



  
AFine=[A{L,1,1} A{L,1,2};
       A{L,2,1} A{L,2,2}];









end



 
    






