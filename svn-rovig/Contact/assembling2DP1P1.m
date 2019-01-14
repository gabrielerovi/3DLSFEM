function A=assembling2DP1P1(mesh,qrule,alpha,beta,coeff1,coeff2,coeff_equilibrium,coeff_constitutive,input_name)

N=mesh.N;
NT=mesh.NT;
node_per_elem=mesh.node_per_elem;

%A=zeros(N,N);
A=sparse(N,N);
if(coeff1==3 && coeff2==3)
for t=1:NT

elem=mesh.elem(t,:);
node=mesh.node(elem,:);

A(elem,elem)=A(elem,elem)+assembling_A33(qrule,node,node_per_elem,coeff_equilibrium,coeff_constitutive,input_name);

end
elseif(coeff1==3 && coeff2==4)
for t=1:NT

elem=mesh.elem(t,:);
node=mesh.node(elem,:);

A(elem,elem)=A(elem,elem)+assembling_A34(qrule,node,node_per_elem,coeff_equilibrium,coeff_constitutive,input_name);

end
elseif(coeff1==4 && coeff2==4)
for t=1:NT

elem=mesh.elem(t,:);
node=mesh.node(elem,:);

A(elem,elem)=A(elem,elem)+assembling_A44(qrule,node,node_per_elem,coeff_equilibrium,coeff_constitutive,input_name);

end

end
