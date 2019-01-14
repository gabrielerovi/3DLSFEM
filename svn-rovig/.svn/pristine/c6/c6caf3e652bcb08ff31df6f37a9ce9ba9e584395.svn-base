
function [Ant,Complementarity,bnt,bnt1]=create_system_LSelasticityContact(parameters,mesh,h,P)

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

%projection
RTCtoRTF= P.RTCtoRTF;
P1CtoP1F= P.P1CtoP1F;



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


A=cell(dim_problem,dim_problem);
A{1,1}=assembling2DRTRT(mesh{L},qrule,alpha,beta,1,1,C_eq,C_const,C_asym,input_name);
A{1,2}=assembling2DRTRT(mesh{L},qrule,alpha,beta,1,2,C_eq,C_const,C_asym,input_name);
A{2,1}=A{1,2}';
A{2,2}=assembling2DRTRT(mesh{L},qrule,alpha,beta,2,2,C_eq,C_const,C_asym,input_name);
A{1,3}=assembling2DRTP1(mesh{L},qrule,alpha,beta,1,3,C_eq,C_const,input_name);
A{1,4}=assembling2DRTP1(mesh{L},qrule,alpha,beta,1,4,C_eq,C_const,input_name);
A{2,3}=assembling2DRTP1(mesh{L},qrule,alpha,beta,2,3,C_eq,C_const,input_name);
A{2,4}=assembling2DRTP1(mesh{L},qrule,alpha,beta,2,4,C_eq,C_const,input_name);
A{3,1}=assembling2DRTP1(mesh{L},qrule,alpha,beta,3,1,C_eq,C_const,input_name);
A{3,2}=assembling2DRTP1(mesh{L},qrule,alpha,beta,3,2,C_eq,C_const,input_name);
A{3,3}=assembling2DP1P1(mesh{L},qrule,alpha,beta,3,3,C_eq,C_const,input_name);
A{3,4}=assembling2DP1P1(mesh{L},qrule,alpha,beta,3,4,C_eq,C_const,input_name);
A{4,1}=assembling2DRTP1(mesh{L},qrule,alpha,beta,4,1,C_eq,C_const,input_name);
A{4,2}=assembling2DRTP1(mesh{L},qrule,alpha,beta,4,2,C_eq,C_const,input_name);
A{4,3}=A{3,4}';
A{4,4}=assembling2DP1P1(mesh{L},qrule,alpha,beta,4,4,C_eq,C_const,input_name);

b1=assemblingb11(mesh{L},qrule,parameters.force1,C_eq,0);
b2=assemblingb11(mesh{L},qrule,parameters.force2,C_eq,0);
b3=assemblingb22(mesh{L},qrule,parameters.fzero,C_eq,C_const);
b4=assemblingb22(mesh{L},qrule,parameters.fzero,C_eq,C_const);

Ant=[];
for mm=1:E_components+N_components
    Araw=[];
    for nn=1:E_components+N_components
        Araw=[Araw, A{mm,nn}];
    end
    Ant=[Ant;Araw];
end
clearvars Araw A

b=[b1;b2;b3;b4];

[M_Normal_Tangent] = MatrixOnGammaCwithNormalTangentComponents(mesh{L},1);

Ant = M_Normal_Tangent * Ant * M_Normal_Tangent';
bnt = M_Normal_Tangent * b;

[Complementarity,bnt_complementarityL] =ComplementarityConditionNT(L,mesh,parameters);

bnt = bnt + bnt_complementarityL;

[Complementarity1,bnt_complementarity1] =ComplementarityConditionNT(1,mesh,parameters);


[M_Normal_Tangent] = MatrixOnGammaCwithNormalTangentComponents(mesh{1},1);

b1=assemblingb11(mesh{1},qrule,parameters.force1,C_eq,0);
b2=assemblingb11(mesh{1},qrule,parameters.force2,C_eq,0);
b3=assemblingb22(mesh{1},qrule,parameters.fzero,C_eq,C_const);
b4=assemblingb22(mesh{1},qrule,parameters.fzero,C_eq,C_const);

b=[b1;b2;b3;b4];

bnt1 = M_Normal_Tangent * b;
bnt1 = bnt1 + h(1)/h(L) * bnt_complementarity1;



clearvars Complementarity1
clearvars M_Normal_Tangent

b_tmp{1}=bnt1;
b_tmp{2}=bnt;
cont=0;
for lev=[1,L]
cont=cont+1;
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
    b_tmp{cont}(jj1)=tmp(1)/coeff;
    b_tmp{cont}(jj2)=tmp(2)/coeff;
end

NcontBC=0;
type_of_dof=1;
for ii=N_remove
    U_xy=[0;0];
    NcontBC=NcontBC+1;
    tmp=add_boundary_bc_elasticity2D(U_xy,type_of_dof,N_label(NcontBC), 0);   
    jj1=ii+2*NE;
    jj2=ii+2*NE+N;
    b_tmp{cont}(jj1)=tmp(1);
    b_tmp{cont}(jj2)=tmp(2);
end
end

bnt=b_tmp{2};
bnt1=b_tmp{1};

end
