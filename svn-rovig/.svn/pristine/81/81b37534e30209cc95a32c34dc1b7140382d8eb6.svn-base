function [MS,MSU,MU,F]=assembling_Contact(qrule,node,face_per_elem,node_per_elem,alpha,beta,Ceq,Cconst,Casym,fx,fy,fz)

[q_point,weights,Volume]=quadrature_points_3D(qrule,node);

number_of_qp=length(q_point(:,1));
[RT0,RT0_div] = phiRT3Dcell(q_point,node);
[P1grad] = P1grad3D(node);

Valss=zeros(number_of_qp,1);
Mss=cell(3,3);
for mm=1:3
for nn=mm:3
Mss{mm,nn}=zeros(face_per_elem,face_per_elem);
Muu{mm,nn}=zeros(node_per_elem,node_per_elem);
end
for nn=1:3
    Msu{mm,nn}=zeros(face_per_elem,face_per_elem);

end
end


 F=zeros(3*face_per_elem,1);

C1=beta*beta;
C2=(2*alpha^2+(alpha+beta)^2);
C3=(alpha^2+2*alpha*(alpha+beta));
C4 = - Cconst*(alpha+beta);
C5 = - 0.5 * Cconst * beta;
C6 = - Cconst * alpha;

   for f1=1:face_per_elem    
        for f2=1:face_per_elem
                        
            %for qp=1:number_of_qp
          Valss(1,1)= sum(weights(:).* ...
                    ((Cconst * C1 +2 * Casym) .* (RT0{f1}(:,2).* RT0{f2}(:,2) + RT0{f1}(:,3).* RT0{f2}(:,3) ) +...
                     Cconst .* C2 .* RT0{f1}(:,1).* RT0{f2}(:,1) ) )+ ... 
                     Ceq .* RT0_div(f1) .* RT0_div(f2) ;
                 
            Valss(1,2,:)= sum(weights(:).* ...
                    (Cconst .* C3 .* RT0{f1}(:,1).* RT0{f2}(:,2) - 2 .* Casym .* RT0{f1}(:,2).* RT0{f2}(:,1)  ));

            Valss(1,3)=sum( weights(:).* ...
                    (Cconst .* C3 .* RT0{f1}(:,1).* RT0{f2}(:,3) - 2 .* Casym .* RT0{f1}(:,3).* RT0{f2}(:,1)  ));
                 
                 
            Valss(2,2)= sum(weights(:).* ...
                    ((Cconst .* C1 +2 .* Casym) .* (RT0{f1}(:,1).* RT0{f2}(:,1) + RT0{f1}(:,3).* RT0{f2}(:,3) ) +...
                     Cconst .* C2 .* RT0{f1}(:,2).* RT0{f2}(:,2) ) ) + ... 
                     Ceq .* RT0_div(f1) .* RT0_div(f2);    
                 
            Valss(2,3)=sum( weights(:).* ...
                    (Cconst .* C3 .* RT0{f1}(:,2).* RT0{f2}(:,3) - 2 .* Casym .* RT0{f1}(:,3).* RT0{f2}(:,2)  ));  
                 
            Valss(3,3)= sum(weights(:).* ...
                    ((Cconst .* C1 +2 .* Casym) .* (RT0{f1}(:,1).* RT0{f2}(:,1) + RT0{f1}(:,2).* RT0{f2}(:,2) ) +...
                     Cconst .* C2 .* RT0{f1}(:,3).* RT0{f2}(:,3)  ) )+ ... 
                     Ceq .* RT0_div(f1) .* RT0_div(f2) ;                 

                 
                 
           Valsu(1,1)= sum(weights(:).* ...
                           (C4 .* RT0{f1}(:,1).*P1grad(f2,1) + C5 .*  ( RT0{f1}(:,2).*P1grad(f2,2) + RT0{f1}(:,3).*P1grad(f2,3)  ) ) );
                 
            Valsu(1,2)=sum( weights(:).* ...
                           (C6 .* RT0{f1}(:,1).*P1grad(f2,2) + C5 .* RT0{f1}(:,2).*P1grad(f2,1) ) );
                   
            Valsu(1,3)= sum(weights(:).* ...
                           (C6 .* RT0{f1}(:,1).*P1grad(f2,3) + C5 .* RT0{f1}(:,3).*P1grad(f2,1) ) );



            Valsu(2,1)=sum( weights(:).* ...
                           (C6 .* RT0{f1}(:,2).*P1grad(f2,1) + C5 .* RT0{f1}(:,1).*P1grad(f2,2) ) );
                 
            Valsu(2,2)= sum(weights(:).* ...
                           (C4 .* RT0{f1}(:,2).*P1grad(f2,2) + C5 .* ( RT0{f1}(:,1).*P1grad(f2,1) + RT0{f1}(:,3).*P1grad(f2,3) ) ) );
                   
            Valsu(2,3)=sum( weights(:).* ...
                           (C6 .* RT0{f1}(:,2).*P1grad(f2,3) + C5 .* RT0{f1}(:,3).*P1grad(f2,2) ) );



            Valsu(3,1)=sum( weights(:).* ...
                           (C6 .* RT0{f1}(:,3).*P1grad(f2,1) + C5 .* RT0{f1}(:,1).*P1grad(f2,3) ) );
                 
            Valsu(3,2)=sum( weights(:).* ...
                           (C6 .* RT0{f1}(:,3).*P1grad(f2,2) + C5 .* RT0{f1}(:,2).*P1grad(f2,3) ) );
                   
            Valsu(3,3)=sum( weights(:).* ...
                           (C4 .* RT0{f1}(:,3).*P1grad(f2,3) + C5 .* ( RT0{f1}(:,1).*P1grad(f2,1) + RT0{f1}(:,2).*P1grad(f2,2) ) ) );
                  
                 
                 

            %end
 
            Valuu(1,1)= Volume * Cconst * (      P1grad(f1,1)*P1grad(f2,1) + 0.5 * P1grad(f1,2)*P1grad(f2,2) + 0.5 * P1grad(f1,3)*P1grad(f2,3));  
                 
            Valuu(1,2)= Volume * Cconst * (0.5 * P1grad(f1,2)*P1grad(f2,1));

            Valuu(1,3)= Volume * Cconst * (0.5 * P1grad(f1,3)*P1grad(f2,1));
                 
            Valuu(2,2)= Volume * Cconst * (0.5 * P1grad(f1,1)*P1grad(f2,1) +     P1grad(f1,2)*P1grad(f2,2) + 0.5 * P1grad(f1,3)*P1grad(f2,3));  
                 
            Valuu(2,3)= Volume * Cconst * (0.5 * P1grad(f1,3)*P1grad(f2,2));

            Valuu(3,3)= Volume * Cconst * (0.5 * P1grad(f1,1)*P1grad(f2,1) + 0.5 * P1grad(f1,2)*P1grad(f2,2) +        P1grad(f1,3)*P1grad(f2,3));  

            
            
          for mm=1:3
              for nn=mm:3
                   Mss{mm,nn}(f1,f2)=Mss{mm,nn}(f1,f2)+Valss(mm,nn)*Volume;  
                   Muu{mm,nn}(f1,f2)=Muu{mm,nn}(f1,f2)+Valuu(mm,nn);                       
              end
              for nn=1:3
                  Msu{mm,nn}(f1,f2)=Msu{mm,nn}(f1,f2)+Valsu(mm,nn)*Volume;
              end
          end
          
        end

        
   end
  
   vec=1:face_per_elem;
   
                F(vec,1)= F(vec,1)+ sum( weights .* (-1.0)*Ceq * RT0_div(vec) .* fx(q_point(:,1),q_point(:,2),q_point(:,3)) ,1)';
                F(vec+face_per_elem,1)= F(vec+face_per_elem,1)+ sum( weights .* (-1.0)*Ceq * RT0_div(vec) .* fy(q_point(:,1),q_point(:,2),q_point(:,3)) ,1)';
                F(vec+2*face_per_elem,1)= F(vec+2*face_per_elem,1)+ sum( weights .* (-1.0)*Ceq * RT0_div(vec) .* fz(q_point(:,1),q_point(:,2),q_point(:,3)) ,1)';
  
                
    F= Volume * F;              

   
          for mm=1:2
              for nn=mm+1:3
                   Mss{nn,mm}=Mss{mm,nn}';
                   Muu{nn,mm}=Muu{mm,nn}';
              end
          end


          
          
          MS=[Mss{1,1} Mss{1,2} Mss{1,3};
              Mss{2,1} Mss{2,2} Mss{2,3};
              Mss{3,1} Mss{3,2} Mss{3,3};];
          
          MU=[Muu{1,1} Muu{1,2} Muu{1,3};
              Muu{2,1} Muu{2,2} Muu{2,3};
              Muu{3,1} Muu{3,2} Muu{3,3};];
          
          MSU=[Msu{1,1} Msu{1,2} Msu{1,3};
               Msu{2,1} Msu{2,2} Msu{2,3};
               Msu{3,1} Msu{3,2} Msu{3,3};];
          
end
