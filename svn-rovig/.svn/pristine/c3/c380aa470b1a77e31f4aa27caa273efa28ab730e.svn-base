function F=assembling_DispVolumeForce(qrule,node,node_per_elem,f)

[q_point,weights,area]=quadrature_points_2D(qrule,node);

number_of_qp=length(q_point(:,1));
P1_basis = phiP12D(q_point,node);
P1grad = P1grad2D(node);
value=zeros(number_of_qp,1);
F=zeros(node_per_elem,1);

   for nn=1:node_per_elem    
            for qp=1:number_of_qp
                value(qp)= weights(qp) * P1_basis(nn,qp) *f(q_point(1),q_point(2));
            end
                 F(nn,1)=F(nn,1)+ area * sum(value);                    
   end
    
end

    