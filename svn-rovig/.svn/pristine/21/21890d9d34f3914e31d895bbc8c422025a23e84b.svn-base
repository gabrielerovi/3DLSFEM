function M=assembling_SigmaU(qrule,node,node_per_elem,face_per_elem,alpha,beta,Ceq,Cconst,Casym)

[q_point,weights,Volume]=quadrature_points_3D(qrule,node);

number_of_qp=length(q_point(:,1));
[RT0,RT0_div] = phiRT3Dcell(q_point,node);
P1grad = P1grad3D(node);
value=zeros(number_of_qp,1);
M=cell(3,3);
for mm=1:3
for nn=1:3
M{mm,nn}=zeros(face_per_elem,face_per_elem);
end
end

C1 = - Cconst*(alpha+beta);
C2 = - 0.5 * Cconst * beta;
C3 = - Cconst * alpha;

   for f1=1:face_per_elem    
        for n2=1:node_per_elem
                        
            for qp=1:number_of_qp

            value(1,1,qp)= weights(qp)* ...
                           (C1 * RT0{f1}(qp,1)*P1grad(n2,1) + C2 *  ( RT0{f1}(qp,2)*P1grad(n2,2) + RT0{f1}(qp,3)*P1grad(n2,3)  ) ) ;
                 
            value(1,2,qp)= weights(qp)* ...
                           (C3 * RT0{f1}(qp,1)*P1grad(n2,2) + C2 * RT0{f1}(qp,2)*P1grad(n2,1) ) ;
                   
            value(1,3,qp)= weights(qp)* ...
                           (C3 * RT0{f1}(qp,1)*P1grad(n2,3) + C2 * RT0{f1}(qp,3)*P1grad(n2,1) ) ;



            value(2,1,qp)= weights(qp)* ...
                           (C3 * RT0{f1}(qp,2)*P1grad(n2,1) + C2 * RT0{f1}(qp,1)*P1grad(n2,2) ) ;
                 
            value(2,2,qp)= weights(qp)* ...
                           (C1 * RT0{f1}(qp,2)*P1grad(n2,2) + C2 * ( RT0{f1}(qp,1)*P1grad(n2,1) + RT0{f1}(qp,3)*P1grad(n2,3) ) ) ;
                   
            value(2,3,qp)= weights(qp)* ...
                           (C3 * RT0{f1}(qp,2)*P1grad(n2,3) + C2 * RT0{f1}(qp,3)*P1grad(n2,2) ) ;



            value(3,1,qp)= weights(qp)* ...
                           (C3 * RT0{f1}(qp,3)*P1grad(n2,1) + C2 * RT0{f1}(qp,1)*P1grad(n2,3) ) ;
                 
            value(3,2,qp)= weights(qp)* ...
                           (C3 * RT0{f1}(qp,3)*P1grad(n2,2) + C2 * RT0{f1}(qp,2)*P1grad(n2,3) ) ;
                   
            value(3,3,qp)= weights(qp)* ...
                           (C1 * RT0{f1}(qp,3)*P1grad(n2,3) + C2 * ( RT0{f1}(qp,1)*P1grad(n2,1) + RT0{f1}(qp,2)*P1grad(n2,2) ) ) ;
                 



            end
          
          for mm=1:3
              for nn=1:3
                   M{mm,nn}(f1,n2)=M{mm,nn}(f1,n2)+sum(value(mm,nn,:))*Volume;  
              end
          end
          
        end
    end
end
