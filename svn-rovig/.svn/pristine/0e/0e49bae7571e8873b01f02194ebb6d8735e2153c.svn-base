function F=assembling_b11(qrule,node,edge_per_elem,f,coeff_equilibrium,coeff_constitutive)

[q_point,weights,area]=quadrature_points_2D(qrule,node);
[RT0,RT0_div] = phiRT2Dcell(q_point,node);
number_of_qp=length(q_point(:,1));
value=zeros(number_of_qp,1);
F=zeros(edge_per_elem,1);

   for nn=1:edge_per_elem    
            for qp=1:number_of_qp
                value(qp)= weights(qp) * (-1.0)*coeff_equilibrium * RT0_div(nn) *f(q_point(qp,1),q_point(qp,2));
            end
                 F(nn,1)=F(nn,1)+ area * sum(value);                    
   end
    
end

    