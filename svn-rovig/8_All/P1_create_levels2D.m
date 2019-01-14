function [A_lev,b_lev,P1CtoP1F]=P1_create_levels2D(A,b,mesh,P1CtoP1F)

L=size(mesh);
L=L(1);
A_lev=cell(L,1);
b_dirichlet=zeros(length(A(:,1)),1);
A_lev{L}=A;
if(L>1)
for lev=L-1:-1:1   
    A_lev{lev}=P1CtoP1F{lev}'*A_lev{lev+1}*P1CtoP1F{lev};
end
end
N_bc=mesh{L}.N_bc;
N_dirichlet=mesh{L}.N_dirichlet;
internal=0;
boundary=0;
type_of_dof=1;
node_remove=[];
for nF=1:mesh{L}.N
    if(N_dirichlet(nF)>0)
        b_dirichlet(nF)=dirichlet_bc(N_bc(nF),type_of_dof); 
        boundary=boundary+1;
        node_remove(boundary)=nF;
    else
% for each internal dof, compute A_{i,:)b_dirichlet(i), that is known and
% can be put in the rhs
        internal=internal+1;
    end
end
rhs=zeros(internal,1);
internal=0;
for nF=1:mesh{L}.N
    if(N_dirichlet(nF)==0)
% for each internal dof, compute A_{i,:)b_dirichlet(i), that is known and
% can be put in the rhs
        internal=internal+1;
        rhs(internal)=A(nF,:)*b_dirichlet;
    end
end

A_lev{L}(node_remove,:)=[];
A_lev{L}(:,node_remove)=[];
b(node_remove)=[];
b_lev=b-rhs;

for lev=1:L-1
    NC_remove=mesh{lev}.N_remove;
    NF_remove=mesh{lev+1}.N_remove;
    
    A_lev{lev}(NC_remove,:)=[];
    A_lev{lev}(:,NC_remove)=[];
    
end

end