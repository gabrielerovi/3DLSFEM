
function [A,b,AFine,P]=create_system_DispElasticity(parameters,mesh,h,P)

% problem coefficients
input_name=parameters.input_name;
qrule=parameters.qrule;
L=length(mesh);
N=mesh{L}.N;
N_components=parameters.N_components;
E_components=parameters.E_components;
N_remove=mesh{L}.N_remove;
N_label=mesh{L}.N_label;
dim_problem=N_components+E_components;

% coefficients
qrule=parameters.qrule;
alpha=parameters.alpha;
beta=parameters.beta;
contact=parameters.contact;

C1=2*parameters.mu;
C2=parameters.lambda;
%projection
P1CtoP1F= P.P1CtoP1F;


if(parameters.h_dependence==false)
    h=ones(L,1);
end

switch N_components
    case 0
        N_remove=[];
        N_label=[];
    otherwise
        N_remove=[mesh{L}.N_remove];
        N_label=mesh{L}.N_label;
end


A=cell(L,dim_problem,dim_problem);

A{L,1,1}=assembling2DP1P1(mesh{L},qrule,alpha,beta,3,3,C1,C2,input_name);
A{L,1,2}=assembling2DP1P1(mesh{L},qrule,alpha,beta,3,4,C1,C2,input_name);
A{L,2,1}=A{L,1,2}';
A{L,2,2}=assembling2DP1P1(mesh{L},qrule,alpha,beta,4,4,C1,C2,input_name);

b1=assemblingDispVolumeForce(mesh{L},qrule,parameters.force1);
b2=assemblingDispVolumeForce(mesh{L},qrule,parameters.force2);

b=[b1;b2];

if(L>1)
for lev=L-1:-1:1   
    
    A{lev,1,1}=P1CtoP1F{lev}'*A{lev+1,1,1}*P1CtoP1F{lev};
    A{lev,1,2}=P1CtoP1F{lev}'*A{lev+1,1,2}*P1CtoP1F{lev};
    A{lev,2,1}=A{lev,1,2}';
    A{lev,2,2}=P1CtoP1F{lev}'*A{lev+1,2,2}*P1CtoP1F{lev};

    
end
end





% add bc to finer level
A{L,1,1}(N_remove,:)=0;
A{L,1,2}(N_remove,:)=0;
A{L,2,1}(N_remove,:)=0;
A{L,2,2}(N_remove,:)=0;


% Dirichlet BC
NcontBC=0;
type_of_dof=1;
for ii=N_remove
    U_xy=[0;0];
    NcontBC=NcontBC+1;
    tmp=add_boundary_bc_elasticity2D(U_xy,type_of_dof,N_label(NcontBC), 0);   
    jj1=ii;
    jj2=ii+N;
    b(jj1)=tmp(1);
    b(jj2)=tmp(2);
    A{L,1,1}(ii,ii)=1;  
    A{L,2,2}(ii,ii)=1; 

end
% Neumann BC
type_of_dof=2;
[dirichlet,n_and_or_t,bool_bc]= boundary_value_bool(type_of_dof);

for ee=1:length(mesh{L}.boundary)
    
    vertices=mesh{L}.boundary(ee,[1,2]);
    node=mesh{L}.node(vertices,[1,2]);
    edgelength=norm(node(1,:)-node(2,:));
    label=mesh{L}.boundary(ee,3);
    
    if(bool_bc(label)==1)
    for nn=1:2   
    if(mesh{L}.N_dirichlet(vertices(nn))==0)    
    b1(vertices(nn))= b1(vertices(nn)) + (0.5 * edgelength) * dirichlet(label,1);   
    b2(vertices(nn))= b2(vertices(nn)) + (0.5 * edgelength) * dirichlet(label,2);
    end
    end
    end

end


% add bc to rougher levels
for lev=1:L-1
N_remove=mesh{lev}.N_remove;
A{lev,1,1}(N_remove,:)=0;
A{lev,1,2}(N_remove,:)=0;
A{lev,2,1}(N_remove,:)=0;
A{lev,2,2}(N_remove,:)=0;

    for ii=N_remove
        A{lev,1,1}(ii,ii)=1;
        A{lev,2,2}(ii,ii)=1;
    end
    
    
end



  
AFine=[A{L,1,1} A{L,1,2} ;
       A{L,2,1} A{L,2,2} ;];




b=b+[b1;b2];




end



 
    






