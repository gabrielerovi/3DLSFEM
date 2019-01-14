function C=assembling2DNDND(mesh,qrule,alpha,beta,coeff1,coeff2,coeff_equilibrium,coeff_constitutive,coeff_symmetry)

N=mesh.N;
NT=mesh.NT;
node_per_elem=mesh.node_per_elem;

%C=zeros(N,N);
C=sparse(N,N);
if(coeff1==1 && coeff2==1)
for t=1:NT

elem=mesh.elem(t,:);
node=mesh.node(elem,:);

C(elem,elem)=C(elem,elem)+assembling_C11(qrule,node,node_per_elem,alpha,beta,coeff_equilibrium,coeff_constitutive,coeff_symmetry);

end
elseif(coeff1==1 && coeff2==2)
    for t=1:NT

elem=mesh.elem(t,:);
node=mesh.node(elem,:);

C(elem,elem)=C(elem,elem)+assembling_C12(qrule,node,node_per_elem,alpha,beta,coeff_equilibrium,coeff_constitutive,coeff_symmetry);

end
% elseif(coeff1==1 && coeff2==3)
% for t=1:NT
% 
% elem=mesh.elem(t,:);
% node=mesh.node(elem,:);
% 
% C(elem,elem)=C(elem,elem)+assembling_C13(qrule,node,node_per_elem,alpha,beta);
% 
% end   
% elseif(coeff1==1 && coeff2==4)
% for t=1:NT
% 
% elem=mesh.elem(t,:);
% node=mesh.node(elem,:);
% 
% C(elem,elem)=C(elem,elem)+assembling_C14(qrule,node,node_per_elem,alpha,beta);
% 
% end
elseif(coeff1==2 && coeff2==2)
for t=1:NT

elem=mesh.elem(t,:);
node=mesh.node(elem,:);

C(elem,elem)=C(elem,elem)+assembling_C22(qrule,node,node_per_elem,alpha,beta,coeff_equilibrium,coeff_constitutive,coeff_symmetry);

end
% elseif(coeff1==2 && coeff2==3)
% for t=1:NT
% 
% elem=mesh.elem(t,:);
% node=mesh.node(elem,:);
% 
% C(elem,elem)=C(elem,elem)+assembling_C23(qrule,node,node_per_elem,alpha,beta);
% 
% end
% elseif(coeff1==2 && coeff2==4)
% for t=1:NT
% 
% elem=mesh.elem(t,:);
% node=mesh.node(elem,:);
% 
% C(elem,elem)=C(elem,elem)+assembling_C24(qrule,node,node_per_elem,alpha,beta);
% 
% end
end
