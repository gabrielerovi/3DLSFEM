function F_loc=rhs_lagrangian(node,qrule)

         [q_point,weights,area]=quadrature_points_2D(qrule,node);
         number_of_qp=length(q_point(:,1));
         P1_basis = phiP12D(q_point,node);
         F_loc=zeros(3,1);
            for k=1:3
                 
                 value=zeros(number_of_qp,1);
                 for qp=1:number_of_qp
                     value(qp)= 1 * weights(qp)*P1_basis(k,qp) ;
                 end
                 q=sum(value)*area;             
                 F_loc(k)=q;            

             end


end
