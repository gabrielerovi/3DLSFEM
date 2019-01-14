function D_value=Assembling_Local_Lumped_Asymmetry_Sigma(qrule,node,edge_per_elem,coeff_symmetry,IsFirstComponentOfStress)

% given an edge e of a triangle T:
% node is the matrix of the coordinates of:
% 1, 2) the vertices of the edge e
% 3) the barycenter of T

[q_point,weights,area]=quadrature_points_2D(qrule,node);

number_of_qp=length(q_point(:,1));
[RT0,RT0_div] = phiRT2Dcell(q_point,node);
value=zeros(number_of_qp,1);

if(IsFirstComponentOfStress)
    %   (sxy-syx) txy
    sign=1;
    tau_component=2;
else
    % -(sxy-syx) tyx
    sign=-1;
    tau_component=1;
end

           for qp=1:number_of_qp
                value(qp)= weights(qp)* coeff_symmetry * sign * ( (RT0{1}(qp,2) - RT0{1}(qp,1) ) * RT0{1}(qp,tau_component)) ;
           end
            
           D_value=sum(value)*area;            

end
