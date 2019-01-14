function [P1_basis] = phiP12D(q_point,node)
% q_point is a matrix m x n
% m = qrule_npoints
% n = dimension of the problem

% node contains the node of the triangle
   
    J(1,1)=node(2,1)-node(1,1);
    J(1,2)=node(3,1)-node(1,1);
    J(2,1)=node(2,2)-node(1,2);
    J(2,2)=node(3,2)-node(1,2);
    detJ=J(1,1)*J(2,2)-J(1,2)*J(2,1);
    
     tmp_matrix=...
    [ node(1,1) node(1,2) 1 ;
      node(2,1) node(2,2) 1 ;
      node(3,1) node(3,2) 1 ;];
 
    f1=[1;0;0]; f2=[0;1;0]; f3=[0;0;1];

    b(1,:)=tmp_matrix\f1 ;
    b(2,:)=tmp_matrix\f2 ;
    b(3,:)=tmp_matrix\f3 ;
    
    
  
   
    qrule_npoints=length(q_point(:,1));
    P1_basis=zeros(3,qrule_npoints);
  
    
    for qp=1:qrule_npoints

        
        % if sign=1/-1, only the shape function c(a,b) has (have) inwardnormal
        for i=1:3
        P1_basis(i,qp)=   b(i,1) * q_point(qp,1)+ b(i,2) * q_point(qp,2) + b(i,3)   ;
        end

    end
    
    
end
        
        