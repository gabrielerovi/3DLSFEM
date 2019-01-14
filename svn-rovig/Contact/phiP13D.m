function [P1_basis,P1grad] = phiP13D(q_point,node)
% q_point is a matrix m x n
% m = qrule_npoints
% n = dimension of the problem

% node contains the node of the triangle
    
     tmp_matrix=...
    [ node(1,1) node(1,2) node(1,3) 1 ;
      node(2,1) node(2,2) node(2,3) 1 ;
      node(3,1) node(3,2) node(3,3) 1 ;
      node(4,1) node(4,2) node(4,3) 1 ;];
 
    f1=[1;0;0;0]; f2=[0;1;0;0]; f3=[0;0;1;0]; f4=[0;0;0;1];

    b(1,:)=tmp_matrix\f1 ;
    b(2,:)=tmp_matrix\f2 ;
    b(3,:)=tmp_matrix\f3 ;
    b(4,:)=tmp_matrix\f4 ;
    
  
   
    qrule_npoints=length(q_point(:,1));
    P1_basis=zeros(4,qrule_npoints);
  
    
    for qp=1:qrule_npoints

        
        for i=1:4
        P1_basis(i,qp)=   b(i,1) * q_point(qp,1)+ b(i,2) * q_point(qp,2) + b(i,3) * q_point(qp,3)   + b(i,4);
        end

    end
    
    
P1grad=[b(1,1:3);b(2,1:3);b(3,1:3);b(4,1:3)];


end
        
        