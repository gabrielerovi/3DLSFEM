function [RT0_basis,RT0_divergence] = phiRT2Dcell(q_point,node)
% q_point is a matrix m x n
% m = qrule_npoints
% n = dimension of the problem

% node contains the node of the triangle
    
    dim=2;

    J(1,1)=node(2,1)-node(1,1);
    J(1,2)=node(3,1)-node(1,1);
    J(2,1)=node(2,2)-node(1,2);
    J(2,2)=node(3,2)-node(1,2);
    detJ=J(1,1)*J(2,2)-J(1,2)*J(2,1);
    
    a(1)=node(1,1); a(2)=node(1,2); 
    b(1)=node(2,1); b(2)=node(2,2); 
    c(1)=node(3,1); c(2)=node(3,2); 
   
    qrule_npoints=length(q_point(:,1));
    Area=0.5 * abs(detJ);
    RT0_basis=zeros(3,qrule_npoints,2);
    if(detJ<0)
        sign=-1;
    else
        sign=+1;
    end
    

    RT0_basis=cell(3,1);
    
    for qp=1:qrule_npoints

        
        % if sign=1/-1, only the shape function c(a,b) has (have) inwardnormal
        for kkk=1:dim
        RT0_basis{1}(qp,kkk)=sign * (+q_point(qp,kkk)-c(kkk))/(2.0*Area) ;
        RT0_basis{2}(qp,kkk)=sign * (+q_point(qp,kkk)-a(kkk))/(2.0*Area) ;
        RT0_basis{3}(qp,kkk)=sign * (-q_point(qp,kkk)+b(kkk))/(2.0*Area) ;
        end

     end
        RT0_divergence(1)= + sign * dim/(2*Area);
        RT0_divergence(2)= + sign * dim/(2*Area);
        RT0_divergence(3)= - sign * dim/(2*Area);
        
end
        
        