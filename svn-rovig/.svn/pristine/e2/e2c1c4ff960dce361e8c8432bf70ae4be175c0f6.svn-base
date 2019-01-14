function M=assembling_A13(qrule,node,edge_per_elem,node_per_elem,alpha,beta,coeff_equilibrium,coeff_constitutive,input_name)

[q_point,weights,area]=quadrature_points_2D(qrule,node);
number_of_qp=length(q_point(:,1));
value=zeros(number_of_qp,1);
M=zeros(edge_per_elem,node_per_elem);

[RT0,RT0_div] = phiRT2Dcell(q_point,node);
P1_basis = phiP12D(q_point,node);
P1grad = P1grad2D(node);

C1=(alpha+beta);
C2=0.5*beta;

   for ee=1:edge_per_elem    
        for nn=1:node_per_elem          
            for qp=1:number_of_qp
                if( strcmp(input_name,'LSelasticity')||strcmp(input_name,'LSelasticityAsymmetric'))
                    value(qp)= weights(qp)*  (-1.0)* ...
                               coeff_constitutive * (C1 * RT0{ee}(qp,1) * P1grad(nn,1) + C2 * RT0{ee}(qp,2)* P1grad(nn,2)) ;
                end
            end
             M(ee,nn)=M(ee,nn)+sum(value)* area;           
        end
   end
    
end

    