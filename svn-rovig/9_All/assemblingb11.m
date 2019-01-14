function b11=assemblingb11(mesh,qrule,f,coeff_equilibrium,coeff_constitutive)

NE=mesh.NE;
NT=mesh.NT;
b11=sparse(NE,1);
edge_per_elem=mesh.edge_per_elem;

for t=1:NT

elem=mesh.elem(t,:);
node=mesh.node(elem,:);
elemE=mesh.elemE(t,:);

b11(elemE)=b11(elemE)+assembling_b11(qrule,node,edge_per_elem,f,coeff_equilibrium,coeff_constitutive);

end


end
