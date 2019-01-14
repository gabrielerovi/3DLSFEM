function grad = P1grad2D  ( node  )
    J(1,1)=node(2,1)-node(1,1);
    J(1,2)=node(3,1)-node(1,1);
    J(2,1)=node(2,2)-node(1,2);
    J(2,2)=node(3,2)-node(1,2);
    Area =0.5*abs(J(1,1)*J(2,2)-J(1,2)*J(2,1) );
    
    tmp_matrix=...
    [1 node(1,1) node(1,2) ;
     1 node(2,1) node(2,2) ;
     1 node(3,1) node(3,2) ;];
 
    f1=[1;0;0]; f2=[0;1;0]; f3=[0;0;1];

    b1=tmp_matrix\f1 ;
    b2=tmp_matrix\f2 ;
    b3=tmp_matrix\f3 ;
    
    grad=[b1(2:3)';b2(2:3)';b3(2:3)'];
    
    end
    