function M=assembling_SigmaSigma(qrule,node,face_per_elem,alpha,beta,Ceq,Cconst,Casym)

[q_point,weights,Volume]=quadrature_points_3D(qrule,node);

number_of_qp=length(q_point(:,1));
[RT0,RT0_div] = phiRT3Dcell(q_point,node);
value=zeros(number_of_qp,1);
M=cell(3,3);
for mm=1:3
for nn=mm:3
M{mm,nn}=zeros(face_per_elem,face_per_elem);
end
end

C1=beta*beta;
C2=(2*alpha^2+(alpha+beta)^2);
C3=(alpha^2+2*alpha*(alpha+beta));
   for f1=1:face_per_elem    
        for f2=1:face_per_elem
                        
            for qp=1:number_of_qp
            value(1,1,qp)= weights(qp)* ...
                    ((Cconst * C1 +2 * Casym) * (RT0{f1}(qp,2)* RT0{f2}(qp,2) + RT0{f1}(qp,3)* RT0{f2}(qp,3) ) +...
                     Cconst * C2 * RT0{f1}(qp,1)* RT0{f2}(qp,1) + ... 
                     Ceq * RT0_div(f1) * RT0_div(f2) ) ;
                 
            value(1,2,qp)= weights(qp)* ...
                    (Cconst * C3 * RT0{f1}(qp,1)* RT0{f2}(qp,2) - 2 * Casym * RT0{f1}(qp,2)* RT0{f2}(qp,1)  );

            value(1,3,qp)= weights(qp)* ...
                    (Cconst * C3 * RT0{f1}(qp,1)* RT0{f2}(qp,3) - 2 * Casym * RT0{f1}(qp,3)* RT0{f2}(qp,1)  );
                 
                 
            value(2,2,qp)= weights(qp)* ...
                    ((Cconst * C1 +2 * Casym) * (RT0{f1}(qp,1)* RT0{f2}(qp,1) + RT0{f1}(qp,3)* RT0{f2}(qp,3) ) +...
                     Cconst * C2 * RT0{f1}(qp,2)* RT0{f2}(qp,2) + ... 
                     Ceq * RT0_div(f1) * RT0_div(f2) ) ;    
                 
            value(2,3,qp)= weights(qp)* ...
                    (Cconst * C3 * RT0{f1}(qp,2)* RT0{f2}(qp,3) - 2 * Casym * RT0{f1}(qp,3)* RT0{f2}(qp,2)  );  
                 
            value(3,3,qp)= weights(qp)* ...
                    ((Cconst * C1 +2 * Casym) * (RT0{f1}(qp,1)* RT0{f2}(qp,1) + RT0{f1}(qp,2)* RT0{f2}(qp,2) ) +...
                     Cconst * C2 * RT0{f1}(qp,3)* RT0{f2}(qp,3) + ... 
                     Ceq * RT0_div(f1) * RT0_div(f2) ) ;                 
                 
            end
          
          for mm=1:3
              for nn=mm:3
                   M{mm,nn}(f1,f2)=M{mm,nn}(f1,f2)+sum(value(mm,nn,:))*Volume;  
              end
          end
          
        end
    end
          for mm=1:2
              for nn=mm+1:3
                   M{nn,mm}=M{mm,nn}';
              end
          end


end