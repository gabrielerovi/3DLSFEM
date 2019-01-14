function M=assembling_A22(qrule,node,edge_per_elem,alpha,beta,coeff_equilibrium,coeff_constitutive,coeff_symmetry,input_name)

[q_point,weights,area]=quadrature_points_2D(qrule,node);

number_of_qp=length(q_point(:,1));
[RT0,RT0_div] = phiRT2Dcell(q_point,node);
value=zeros(number_of_qp,1);
M=zeros(edge_per_elem,edge_per_elem);

C1=beta*beta;
C2=(alpha^2+(alpha+beta)^2);

   for e1=1:edge_per_elem    
        for e2=1:edge_per_elem
                        
            for qp=1:number_of_qp
                if( strcmp(input_name,'LSstressblock') || strcmp(input_name,'LSelasticity')||strcmp(input_name,'LSelasticityAsymmetric'))
                value(qp)= weights(qp)* ...
                    (coeff_symmetry * RT0{e1}(qp,1) * RT0{e2}(qp,1) + ...
                     coeff_constitutive * C1 * RT0{e1}(qp,1)* RT0{e2}(qp,1) + ...
                     coeff_constitutive * C2 * RT0{e1}(qp,2)* RT0{e2}(qp,2) + ... 
                     coeff_equilibrium * RT0_div(e1) * RT0_div(e2) ) ;
                end
            end
            
                 coeff_symmetry_average=0;
                 intxye1=0;
                 intxye2=0;
                 for qp=1:number_of_qp
                     intxye1=weights(qp)*RT0{e1}(qp,1);
                     intxye2=weights(qp)*RT0{e2}(qp,1);
                 end
                 
                 
                 M(e1,e2)=M(e1,e2)+sum(value)*area + coeff_symmetry_average * intxye1*intxye2;       
        end
    end
end
