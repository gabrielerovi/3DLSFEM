function b22=assemblingb22(mesh,qrule,f,coeff_equilibrium,coeff_constitutive)

N=mesh.N;
NT=mesh.NT;
%b22=zeros(N,1);
b22=sparse(N,1);
node_per_elem=mesh.node_per_elem;

for t=1:NT

elem=mesh.elem(t,:);
node=mesh.node(elem,:);

b22(elem)=b22(elem)+assembling_b22(qrule,node,node_per_elem,f,coeff_equilibrium,coeff_constitutive);

end


end
