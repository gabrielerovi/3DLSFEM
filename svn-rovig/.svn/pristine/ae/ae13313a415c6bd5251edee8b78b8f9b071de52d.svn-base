function [M_Normal_Tangent] = MatrixOnGammaCwithNormalTangentComponentsLOCAL(mesh)
% We want to express each variable on GammaC wrt its normal and tangent
% components, i.e; x= x_n + x_t
% EDGE variables:
% sigma n = s_n n + s_t t, i.e 
% c_1=(s_n n1  + s_t n2) / (phidotn)
% c_2=(s_n n2  - s_t n1) / (phidotn)
% NODE variables:
% u n = u_n n + u_t t, i.e 
% c_1=(u_n n1  + u_t n2)
% c_2=(u_n n2  - u_t n1)
% So basically we are writing, 


N=mesh.N;
NF=mesh.NF;
node=mesh.node;
face=mesh.face;
totdofs=3 * (NF + N);


M_Normal_Tangent=speye(totdofs,totdofs);
if(~isempty(mesh.F_contact))
for bb=mesh.F_contact
        
    normal_F=normal_contact(mean(node(face(bb,:),:)),parameters);
    phidotn=phi_dot_n(mesh,bb);
        
    [xn,H]=HouseHolderTransformation(sparse(3,3), normal_F);  
    M_Normal_Tangent( [bb, bb+NF, bb+2*NF] , :)=0;
    M_Normal_Tangent( [bb, bb+NF, bb+2*NF] , [bb, bb+NF, bb+2*NF] )= H*sign(phidotn);
        
end
end

if(~isempty(mesh.N_contact))
for mm=1:length(mesh.N_contact)

        nn=mesh.N_contact(mm);
        normalN=normal_contact(node(nn,:),parameters);
        
            if(mesh.N_dirichlet(nn)==0)
            for jj=1:3
                 nnode(jj) = 3 * NF + (jj-1)*N + nn;
            end
            [xn,H]=HouseHolderTransformation(sparse(3,3), normalN);
            
            M_Normal_Tangent( [nnode(1), nnode(2), nnode(3)] , : )=0;
            M_Normal_Tangent( [nnode(1), nnode(2), nnode(3)] , [nnode(1), nnode(2), nnode(3)] )=H;    
            end
end
end

end















