function [M_Normal_Tangent] = MatrixOnGammaCwithNormalTangentComponents(mesh,householder)
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
NE=mesh.NE;
E_bc=mesh.E_bc;
edge=mesh.edge;
node=mesh.node;
% we consider normal_node_contact instead of normal_node
% because it is better in case of angles
% normal_node_contact=mesh.normal_node_contact;

% normal_edge=mesh.normal_edge_contact;

[dirichlet,n_and_or_t,bool_bc]= boundary_value_bool(1);


[dirichlet_E,n_and_or_t,bool_bc]= boundary_value_bool(2);

M_Normal_Tangent=speye(2 * NE + 2 * N, 2 * NE + 2 * N);
M_Normal_TangentT=speye(2 * NE + 2 * N, 2 * NE + 2 * N);

for bb=mesh.E_contact
    

        for nn=1:2
        nodeN(nn)=edge(bb,nn);
%         normalN(:,nn)=normal_node_contact{edge(bb,nn)};
        normalN(:,nn)=normal_contact(node(edge(bb,nn),:),parameters);
        end
        %normalE=normal_edge{bb};
        normal_E=normal_contact(mean(node(edge(bb,:),:)),parameters);
        
        
        phidotn=phi_dot_n(mesh,bb);
        
         M_Normal_Tangent( [bb, bb+NE] , :)=0;
         M_Normal_Tangent( [bb, bb+NE] , [bb, bb+NE] )= 1.0/phidotn * [ normal_E(1)  normal_E(2);
                                                                       normal_E(2) -normal_E(1);];
        if(householder)
        [xn,H]=HouseHolderTransformation(zeros(2,2), normal_E);  
        M_Normal_Tangent( [bb, bb+NE] , :)=0;
        M_Normal_Tangent( [bb, bb+NE] , [bb, bb+NE] )= H*sign(phidotn);
        M_Normal_TangentT( [bb, bb+NE] , :)=0;
        M_Normal_TangentT( [bb, bb+NE] , [bb, bb+NE] )= H'*sign(phidotn);
        end
        
        for nn=1:2
            nn1 = 2 * NE + nodeN(nn);
            nn2 = 2 * NE + N + nodeN(nn);
        M_Normal_Tangent( [nn1, nn2] , : )=0;
        M_Normal_Tangent( [nn1, nn2] , [nn1, nn2] )=  [ normalN(1,nn)  normalN(2,nn);
                                                        normalN(2,nn) -normalN(1,nn)];
        if(householder)
        [xn,H]=HouseHolderTransformation(zeros(2,2), normalN(:,nn));
        M_Normal_Tangent( [nn1, nn2] , : )=0;
        M_Normal_Tangent( [nn1, nn2] , [nn1, nn2] )=H;
        M_Normal_TangentT( [nn1, nn2] , : )=0;
        M_Normal_TangentT( [nn1, nn2] , [nn1, nn2] )=H';
        end
        
        end   
        
        
        
                   
end





end

















