
function [Ant,Complementarity,bnt,bnt1]=create_system_LSelasticityContact(parameters,mesh,h,P)

% problem coefficients
qrule=parameters.qrule;
L=length(mesh);
NF=mesh{L}.NF;
N=mesh{L}.N;
NT=mesh{L}.NT;
dim=parameters.dim;
face_per_elem=mesh{L}.face_per_elem;
node_per_elem=mesh{L}.node_per_elem;

% coefficients 
Ceq=parameters.C_eq;
Cconst=parameters.C_const;
Casym=parameters.C_asym;
qrule=parameters.qrule;
alpha=parameters.alpha;
beta=parameters.beta;






Ant=sparse(dim*(N+NF),dim*(N+NF));
bnt=sparse(dim*(N+NF),1);

C(1)=beta*beta;
C(2)=Cconst*(2*alpha^2+(alpha+beta)^2);
C(3)=(alpha^2+2*alpha*(alpha+beta));
C(4) = - Cconst*(alpha+beta);
C(5) = - 0.5 * Cconst * beta;
C(6) = - Cconst * alpha;
C(7)=   Cconst;
C(8)=(Cconst * C(3) - 2 * Casym);
C(9)=(Cconst * C(1) +2 * Casym);
C(10)=Cconst*C(3);
C(11)=-2*Casym;
  [RT0_mass,RT0_div,P1_Grad,RT0_Grad,RT0_divergence,signRT0_divergence,q_point,weights]=assemblingreference();
  



for t=1:NT 
    elem=full(mesh{L}.elem(t,:));
    node=mesh{L}.node(elem,:);
    elemF=full(mesh{L}.elemF(t,:));
    [MS,MSU,MU]=assembling_Contact3(C,RT0_mass,RT0_div,P1_Grad,RT0_Grad,RT0_divergence,node,face_per_elem,Ceq);
    vecN=[elem,elem+N,elem+2*N]+3*NF;
    vecF=[elemF,elemF+NF,elemF+2*NF];
    Ant(vecF,vecF) =    Ant(vecF ,vecF ) +    MS;
    Ant(vecN,vecN) =    Ant(vecN ,vecN)  +    MU;
    Ant(vecF,vecN) =    Ant(vecF,vecN)   +    MSU;
    Ant(vecN,vecF) =    Ant(vecN,vecF)   +    MSU';   
end




if(parameters.istheexternalforcenonzero)
for t=1:NT 
    elem=full(mesh{L}.elem(t,:));
    node=mesh{L}.node(elem,:);
    elemF=full(mesh{L}.elemF(t,:));
    [bss]=assembling_ContactRhs(RT0_divergence,node,face_per_elem,Ceq,parameters.force1,parameters.force2,parameters.force3,q_point,weights,parameters.istheexternalforcenonzero);
    vecF=[elemF,elemF+NF,elemF+2*NF];
    bnt(vecF) = bnt(vecF) + bss;
    
end
end



















[M_Normal_Tangent] = MatrixOnGammaCwithNormalTangentComponents(mesh{L});

Ant = M_Normal_Tangent * Ant * M_Normal_Tangent';
bnt = M_Normal_Tangent * bnt;

[Complementarity,bnt_complementarityL] =ComplementarityConditionNT3D(L,mesh,parameters);

bnt = bnt + bnt_complementarityL;
% we enforce frictionless conditions
zerofriction=[mesh{L}.F_contact+mesh{L}.NF,mesh{L}.F_contact+mesh{L}.NF*2];
bnt(zerofriction)=0;







bnt1=sparse(dim*(mesh{1}.N+mesh{1}.NF),1);



if(parameters.istheexternalforcenonzero)
for t=1:mesh{1}.NT
    elem=full(mesh{1}.elem(t,:));
    node=mesh{1}.node(elem,:);
    elemF=full(mesh{1}.elemF(t,:));
    [bss]=assembling_ContactRhs(RT0_divergence,node,face_per_elem,Ceq,parameters.force1,parameters.force2,parameters.force3,q_point,weights,parameters.istheexternalforcenonzero);
    vecF=[elemF,elemF+NF,elemF+2*NF];
    bnt1(vecF) = bnt1(vecF) + bss;
    
end
end












[Complementarity1,bnt_complementarity1] =ComplementarityConditionNT3D(1,mesh,parameters);


[M_Normal_Tangent] = MatrixOnGammaCwithNormalTangentComponents(mesh{1});


bnt1 = M_Normal_Tangent * bnt1;
bnt1 = bnt1 + h(1)/h(L) * bnt_complementarity1;

% we enforce frictionless conditions
zerofriction=[mesh{1}.F_contact+mesh{1}.NF,mesh{1}.F_contact+mesh{1}.NF*2];
bnt1(zerofriction)=0;

clearvars Complementarity1
clearvars M_Normal_Tangent

b_tmp{1}=bnt1;
b_tmp{2}=bnt;
cont=0;
for lev=[1,L]
cont=cont+1;
F_remove=mesh{lev}.F_remove;
N_remove=mesh{lev}.N_remove;
F_label=mesh{lev}.F_label;
N_label=mesh{lev}.N_label;
NF=mesh{lev}.NF;
N=mesh{lev}.N;

FcontBC=0;
type_of_dof=3;
for ii=F_remove
    FcontBC=FcontBC+1;
    tmp=add_boundary_bc_elasticity3D(type_of_dof,F_label(FcontBC),dim);  
    coeff = phi_dot_n(mesh{lev},ii);
    jj1=ii;
    jj2=ii+NF;    
    jj3=ii+2*NF;    
    b_tmp{cont}(jj1)=tmp(1)/coeff;
    b_tmp{cont}(jj2)=tmp(2)/coeff;
    b_tmp{cont}(jj3)=tmp(3)/coeff;
end

NcontBC=0;
type_of_dof=1;
for ii=N_remove
    NcontBC=NcontBC+1;
    tmp=add_boundary_bc_elasticity3D(type_of_dof,N_label(NcontBC),dim);   
    jj1=ii+3*NF;
    jj2=ii+3*NF+N;
    jj3=ii+3*NF+2*N;
    b_tmp{cont}(jj1)=tmp(1);
    b_tmp{cont}(jj2)=tmp(2);
    b_tmp{cont}(jj3)=tmp(3);
end
end

bnt=b_tmp{2};
bnt1=b_tmp{1};

end
