function F_loc=rhs_raviart(node,qrule,f1,f2)

         [q_point,weights,area]=quadrature_points_2D(qrule,node);
         number_of_qp=length(q_point(:,1));
         [RT0_basis,RT0_divergence] = phiRT2D(q_point,node)    ;     
          F_loc=zeros(3,1);
            for k=1:3
                 
                 value=zeros(number_of_qp,1);
                 for qp=1:number_of_qp
                     value(qp)= weights(qp) * ...
                     ( RT0_basis(k,qp,1) * f1(q_point(1),q_point(2))+ RT0_basis(k,qp,2) * f2(q_point(1),q_point(2)) );
                 end
                 q=sum(value)*area;             
                 F_loc(k)=q;            

             end


end
