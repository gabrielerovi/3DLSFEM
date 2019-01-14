function M = local_div_matrix ( node  )

    J(1,1)=node(2,1)-node(1,1);
    J(1,2)=node(3,1)-node(1,1);
    J(2,1)=node(2,2)-node(1,2);
    J(2,2)=node(3,2)-node(1,2);
    detJ=J(1,1)*J(2,2)-J(1,2)*J(2,1);
    abs_detJ=abs(detJ);
    Area =0.5 * abs_detJ;

    M=zeros(3,3);
    
    Div=2.0/abs_detJ *[1;1;-1];
    
    M=Div*Div';
    M=M*Area;
    end
    