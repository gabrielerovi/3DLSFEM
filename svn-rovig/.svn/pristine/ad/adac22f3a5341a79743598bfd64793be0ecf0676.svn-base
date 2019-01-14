function b11=assemblingDispVolumeForce(mesh,qrule,f,coeff_equilibrium,coeff_constitutive)
N=mesh.N;
NT=mesh.NT;
b11=sparse(N,1);
node_per_elem=mesh.node_per_elem;

for t=1:NT

elem=mesh.elem(t,:);
node=mesh.node(elem,:);

b11(elem)=b11(elem)+assembling_DispVolumeForce(qrule,node,node_per_elem,f);

end


end
