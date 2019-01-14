
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
  




  
% for t=1:NT 
%     elem=full(mesh{L}.elem(t,:));
%     node=mesh{L}.node(elem,:);
%     elemF=full(mesh{L}.elemF(t,:));
%     [MS,MSU,MU]=assembling_Contact3(C,RT0_mass,RT0_div,P1_Grad,RT0_Grad,RT0_divergence,node,face_per_elem,Ceq);
%     vecN=[elem,elem+N,elem+2*N]+3*NF;
%     vecF=[elemF,elemF+NF,elemF+2*NF];
%     Ant(vecF,vecF) =    Ant(vecF ,vecF ) +    MS;
%     Ant(vecN,vecN) =    Ant(vecN ,vecN)  +    MU;
%     Ant(vecF,vecN) =    Ant(vecF,vecN)   +    MSU;
%     Ant(vecN,vecF) =    Ant(vecN,vecF)   +    MSU';   
%     
%     
% end



  Val=zeros(24*24,NT);
  Jtot=zeros(24*24,NT);
  Itot=zeros(24*24,NT);
for t=1:NT 
    elem=full(mesh{L}.elem(t,:));
    node=mesh{L}.node(elem,:);
    elemF=full(mesh{L}.elemF(t,:));
 

    Val(:,t)=assembling_Contact2(C,RT0_mass,RT0_div,P1_Grad,RT0_Grad,RT0_divergence,node,face_per_elem,Ceq);
    J=[elemF, elemF+NF, elemF+2*NF, elem+3*NF, elem+3*NF+N, elem+3*NF+2*N];
    I=J';
    J1=repmat(J,length(J),1);
    Jtot(:,t)=J1(:);
    I1=repmat(I,1,length(I));
    Itot(:,t)=I1(:);
    
    
end  
  
Itot=Itot(:);
Jtot=Jtot(:);
Val=Val(:);
Lval=3*(N+NF);
Ant=sparse(Itot,Jtot,Val,Lval,Lval);


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


% ww=unique([mesh{L}.F_contact+mesh{L}.NF,mesh{L}.F_contact+2*mesh{L}.NF, find(mesh{L}.F_dirichlet),(find(mesh{L}.N_dirichlet)+mesh{L}.NF*3)]);
% Ant(ww,:)=0;
% Ant(:,ww)=0;
% 
% for ii=1:length(ww)
%    
% Ant(ww(ii),ww(ii))=1;
% end
% eigsA=eig(full(Ant));
% eigsC=eig(full(Complementarity));
% eigsAC=eig(full(Ant+Complementarity));
% min(eigsA)
% min(eigsC)
% min(eigsAC)
% 
% 
% 
% 
% 
% 
% F_contact_bool=zeros(NF,1);
% N_contact_bool=zeros(N,1);
% F_contact_bool(mesh{L}.F_contact)=1;
% N_contact_bool(mesh{L}.N_contact)=1;
% Antprova=sparse(dim*(N+NF),dim*(N+NF));
% Cconstitutive(1)=beta*beta;
% Cconstitutive(2)=Cconst*(2*alpha^2+(alpha+beta)^2);
% Cconstitutive(3)=(alpha^2+2*alpha*(alpha+beta));
% Cconstitutive(4) = - Cconst*(alpha+beta);
% Cconstitutive(5) = - 0.5 * Cconst * beta;
% Cconstitutive(6) = - Cconst * alpha;
% Cconstitutive(7)=   Cconst;
% Cconstitutive(8)=(Cconst * C(3) );
% Cconstitutive(9)=(Cconst * C(1));
% Cconstitutive(10)=Cconst*C(3);
% Cconstitutive(11)=0;
% 
% for t=1:NT 
%     elem=full(mesh{L}.elem(t,:));
%     node=mesh{L}.node(elem,:);
%     
% %     ref=0.0001;
% %     node=[0 0 0; ref  0 0; 0 ref 0; 0 0 ref];
%     elemF=full(mesh{L}.elemF(t,:));
%     [MS,MSU,MU]=assembling_Contact3(Cconstitutive,RT0_mass,RT0_div,P1_Grad,RT0_Grad,RT0_divergence,node,face_per_elem,0);
%     vecN=[elem,elem+N,elem+2*N]+3*NF;
%     vecF=[elemF,elemF+NF,elemF+2*NF];
%  
%     finalvalue=0;
% M_Normal_Tangent=speye(24);    
% grid.F_contact=[];
% grid.N_contact=[];
% for ii=1:4
%     if(F_contact_bool(elemF(ii)))
%      grid.F_contact=[grid.F_contact;ii];   
%     end
%     if(N_contact_bool(elem(ii)))
%         grid.N_contact=[grid.N_contact;ii];
%     end
% end
% 
% if(~isempty(grid.F_contact) || ~isempty(grid.N_contact))
% for ii=1:4
% grid.normal_face{ii}=mesh{L}.normal_face{elemF(ii)};
% grid.N_dirichlet(ii)=mesh{L}.N_dirichlet(elem(ii));
% end
% grid.N=4;
% grid.NF=4;
% grid.node=node;
% grid.face=[2 3 4;1 3 4; 1 2 4; 1 2 3];
% 
% [M_Normal_Tangent] = MatrixOnGammaCwithNormalTangentComponentsLOCAL(grid);
% end
% 
%   vec=[vecF,vecN];
%   M=[MS,MSU;MSU',MU];
%   M=M_Normal_Tangent*M*M_Normal_Tangent';
%   
%   if(~isempty(grid.F_contact)&&~isempty(grid.N_contact)&&length(grid.F_contact)<=1 && length(grid.N_contact)<=3)
%       vecloc=[grid.F_contact;grid.NF*3+grid.N_contact];
%       Mnt=M(vecloc,vecloc);
%       mineigM=min(eig(Mnt));
%       maxeigM=max(eig(Mnt));
%       nodes_coord_of_face=grid.node;
%       nodes_coord_of_face(grid.F_contact,:)=[];
%       A=AreaTriangle(nodes_coord_of_face) ;
% 
%       finalvalue= parameters.C_contact* A*0.5774/mineigM;
%       finalvalue= parameters.C_contact* A*0.5774/maxeigM;
%       M(vecloc,vecloc)*finalvalue
%       finalvalue
%       if(mineigM<0.00001)
%           problemaaaaaaaaa=1
%       end
%   elseif(length(grid.F_contact)>1 ||  length(grid.N_contact)==4)
%       problemaaaaaaaaa=2
%   end
%   
% 
% 
%     [MS,MSU,MU]=assembling_Contact3(C,RT0_mass,RT0_div,P1_Grad,RT0_Grad,RT0_divergence,node,face_per_elem,Ceq);
%   M=[MS,MSU;MSU',MU];
%   
%   M=M_Normal_Tangent*M*M_Normal_Tangent';
%   
% % 
% %       if(finalvalue>0)
% %       Z=sparse(24,24);
% %       Z(vecloc,vecloc)=M(vecloc,vecloc);
% %       Antprova(vec,vec)=Antprova(vec,vec)+M;
% %       end
% 
%    M(vecloc,vecloc)=M(vecloc,vecloc)*(1+finalvalue);
%    Antprova(vec,vec)=Antprova(vec,vec)+M;
% end
% 
% 
% 
% % eigsAnt=eig(full(Ant+Complementarity));
% % eigsAnt(find(eigsAnt<0))
% % eigsAntprova=eig(full(Antprova+Complementarity));
% % eigsAntprova(find(eigsAntprova<0))
% % 
% % min(real(eigsAnt))
% % min(real(eigsAntprova))
% 
% 
% Ant=Antprova;



















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
