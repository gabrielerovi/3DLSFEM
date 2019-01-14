function M=assembling_C12(qrule,node,node_per_elem,alpha,beta,coeff_equilibrium,coeff_constitutive,coeff_symmetry)

[q_point,weights,area]=quadrature_points_2D(qrule,node);

number_of_qp=length(q_point(:,1));
P1_basis = phiP12D(q_point,node);
P1grad = P1grad2D(node);
value=zeros(number_of_qp,1);
M=zeros(node_per_elem,node_per_elem);
C1=2*coeff_symmetry;
C2=-2*coeff_constitutive*alpha*(alpha+beta);

   for e1=1:node_per_elem    
        for e2=1:node_per_elem
                        
            for qp=1:number_of_qp
                value(qp)= weights(qp)*(C1 * P1grad(e1,1)*P1grad(e2,2) + C2*P1grad(e1,2)*P1grad(e2,1));
            end
                 M(e1,e2)=M(e1,e2)+sum(value)*area;            
        end
   end
   
    
end

    