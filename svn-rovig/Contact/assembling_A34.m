function M=assembling_A34(qrule,node,node_per_elem,C1,C2,input_name)

[q_point,weights,area]=quadrature_points_2D(qrule,node);

number_of_qp=length(q_point(:,1));
P1_basis = phiP12D(q_point,node);
P1grad = P1grad2D(node);
value=zeros(number_of_qp,1);
M=zeros(node_per_elem,node_per_elem);

   for e1=1:node_per_elem    
        for e2=1:node_per_elem
%             for qp=1:number_of_qp
%                 % for now = 0
%                 value(qp)= 0;
%             end
            if( strcmp(input_name,'LSelasticity')||strcmp(input_name,'LSelasticityAsymmetric') )
                    q= area * C2 * 0.5 * P1grad(e2,1)*P1grad(e1,2);
            elseif( strcmp(input_name,'DispElasticity'))
                    q= area * (C1 * 0.5 * P1grad(e2,1)*P1grad(e1,2) +...
                               C2 * P1grad(e1,1)*P1grad(e2,2)) ;
            end
            %q=q * area;
                     M(e1,e2)=M(e1,e2)+q;    
            

        end
   end
    
end

    