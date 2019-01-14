function A=assembling2DRTRT(mesh,qrule,alpha,beta,coeff1,coeff2,coeff_equilibrium,coeff_constitutive,coeff_symmetry,input_name)

N=mesh.N;
NT=mesh.NT;
NE=mesh.NE;
%A=zeros(NE,NE);
A=sparse(NE,NE);
edge_per_elem=mesh.edge_per_elem;

if(coeff1==1 && coeff2==1)
for t=1:NT

elem=mesh.elem(t,:);
node=mesh.node(elem,:);
elemE=mesh.elemE(t,:);

A(elemE,elemE)=A(elemE,elemE)+assembling_A11(qrule,node,edge_per_elem,alpha,beta,coeff_equilibrium,coeff_constitutive,coeff_symmetry,input_name);

end
elseif(coeff1==1 && coeff2==2)
 for t=1:NT

elem=mesh.elem(t,:);
node=mesh.node(elem,:);
elemE=mesh.elemE(t,:);

A(elemE,elemE)=A(elemE,elemE)+assembling_A12(qrule,node,edge_per_elem,alpha,beta,coeff_equilibrium,coeff_constitutive,coeff_symmetry,input_name);

end   
 elseif(coeff1==2 && coeff2==2)
 for t=1:NT

elem=mesh.elem(t,:);
node=mesh.node(elem,:);
elemE=mesh.elemE(t,:);

A(elemE,elemE)=A(elemE,elemE)+assembling_A22(qrule,node,edge_per_elem,alpha,beta,coeff_equilibrium,coeff_constitutive,coeff_symmetry,input_name);

 end   

 
 

end
