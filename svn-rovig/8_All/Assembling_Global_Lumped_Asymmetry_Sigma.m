function [DiagAsymmetry1,DiagAsymmetry2]=Assembling_Global_Lumped_Asymmetry_Sigma(mesh,qrule,coeff_symmetry)

N=mesh.N;
NT=mesh.NT;
NE=mesh.NE;
%A=zeros(NE,NE);
DiagAsymmetry1=zeros(NE,1);
DiagAsymmetry2=zeros(NE,1);

edge_per_elem=mesh.edge_per_elem;

reordernode=[3 1 2];
for t=1:NT

elem=mesh.elem(t,:);
node=mesh.node(elem,:);
elemE=mesh.elemE(t,:);

barycenter=mean(node);

for ee=1:edge_per_elem
    
nodeLumped=node;
nodeLumped(reordernode(ee),:)=[];
nodeLumped(end+1,:)=barycenter;
DiagAsymmetry1(elemE(ee))=DiagAsymmetry1(elemE(ee))+Assembling_Local_Lumped_Asymmetry_Sigma(qrule,nodeLumped,edge_per_elem,coeff_symmetry,1);
DiagAsymmetry2(elemE(ee))=DiagAsymmetry2(elemE(ee))+Assembling_Local_Lumped_Asymmetry_Sigma(qrule,nodeLumped,edge_per_elem,coeff_symmetry,0);
end


end   

 
 

end
