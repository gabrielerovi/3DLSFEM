function [P1grad] = P1grad3D(node)
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
    
  
     
P1grad=[b(1,1:3);b(2,1:3);b(3,1:3);b(4,1:3)];



end
        
        