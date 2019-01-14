function M = local_mass_matrix ( node  )

    J(1,1)=node(2,1)-node(1,1);
    J(1,2)=node(3,1)-node(1,1);
    J(2,1)=node(2,2)-node(1,2);
    J(2,2)=node(3,2)-node(1,2);
    
    M=zeros(3,3);
    detJ=J(1,1)*J(2,2)-J(1,2)*J(2,1);
    a1=1.0/(12.0*abs(detJ));
    a2=1.0/(4.0*abs(detJ));
    cxx=J(1,1)^2+J(2,1)^2;
    cyy=J(1,2)^2+J(2,2)^2;
    cxy=J(1,2)*J(1,1)+J(2,1)*J(2,2);
    
    
    M(1,2)=cxx * (a1)  + cyy * (-a1) + cxy * (-a1);
    M(1,3)=cxx * (a1)  + cyy * (a1) + cxy * (-a2);
    M(2,3)=cxx * (a1)  + cyy * (-a1) + cxy * (a1);
    M=M+M';

    M(1,1)=cxx * (a1)  + cyy * (a2) + cxy * (-a2);
    M(2,2)=cxx * (a1)  + cyy * (a1) + cxy * (a1);
    M(3,3)=cxx * (a2)  + cyy * (a1) + cxy * (-a2);

    end
    