function A=assembling2DRTP1(mesh,qrule,alpha,beta,coeff1,coeff2,coeff_equilibrium,coeff_constitutive,input_name)

N=mesh.N;
NT=mesh.NT;
NE=mesh.NE;
%A=zeros(NE,N);
if(coeff1<coeff2)
A=sparse(NE,N);
else
A=sparse(N,NE);
end
edge_per_elem=mesh.edge_per_elem;
node_per_elem=mesh.node_per_elem;

if(coeff1==1 && coeff2==3)
for t=1:NT

elem=mesh.elem(t,:);
node=mesh.node(elem,:);
elemE=mesh.elemE(t,:);

A(elemE,elem)=A(elemE,elem)+assembling_A13(qrule,node,edge_per_elem,node_per_elem,alpha,beta,coeff_equilibrium,coeff_constitutive,input_name);

end
elseif(coeff1==1 && coeff2==4)
for t=1:NT

elem=mesh.elem(t,:);
node=mesh.node(elem,:);
elemE=mesh.elemE(t,:);

A(elemE,elem)=A(elemE,elem)+assembling_A14(qrule,node,edge_per_elem,node_per_elem,alpha,beta,coeff_equilibrium,coeff_constitutive,input_name);

end

elseif(coeff1==2 && coeff2==3)
for t=1:NT

elem=mesh.elem(t,:);
node=mesh.node(elem,:);
elemE=mesh.elemE(t,:);

A(elemE,elem)=A(elemE,elem)+assembling_A23(qrule,node,edge_per_elem,node_per_elem,alpha,beta,coeff_equilibrium,coeff_constitutive,input_name);

end
elseif(coeff1==2 && coeff2==4)
for t=1:NT

elem=mesh.elem(t,:);
node=mesh.node(elem,:);
elemE=mesh.elemE(t,:);

A(elemE,elem)=A(elemE,elem)+assembling_A24(qrule,node,edge_per_elem,node_per_elem,alpha,beta,coeff_equilibrium,coeff_constitutive,input_name);

end




elseif(coeff1==3 && coeff2==1)
for t=1:NT

elem=mesh.elem(t,:);
node=mesh.node(elem,:);
elemE=mesh.elemE(t,:);

A(elem,elemE)=A(elem,elemE)+assembling_A31(qrule,node,edge_per_elem,node_per_elem,alpha,beta,coeff_equilibrium,coeff_constitutive,input_name);

end
elseif(coeff1==3 && coeff2==2)
for t=1:NT

elem=mesh.elem(t,:);
node=mesh.node(elem,:);
elemE=mesh.elemE(t,:);

A(elem,elemE)=A(elem,elemE)+assembling_A32(qrule,node,edge_per_elem,node_per_elem,alpha,beta,coeff_equilibrium,coeff_constitutive,input_name);

end
elseif(coeff1==4 && coeff2==1)
for t=1:NT

elem=mesh.elem(t,:);
node=mesh.node(elem,:);
elemE=mesh.elemE(t,:);

A(elem,elemE)=A(elem,elemE)+assembling_A41(qrule,node,edge_per_elem,node_per_elem,alpha,beta,coeff_equilibrium,coeff_constitutive,input_name);

end
elseif(coeff1==4 && coeff2==2)
for t=1:NT

elem=mesh.elem(t,:);
node=mesh.node(elem,:);
elemE=mesh.elemE(t,:);

A(elem,elemE)=A(elem,elemE)+assembling_A42(qrule,node,edge_per_elem,node_per_elem,alpha,beta,coeff_equilibrium,coeff_constitutive,input_name);

end

end
