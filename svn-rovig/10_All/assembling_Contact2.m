function [MS,MSU,MU,F]=assembling_Contact2(qrule,node,face_per_elem,node_per_elem,alpha,beta,Ceq,Cconst,Casym,fx,fy,fz)


xw = TriGaussPoints3D(qrule);
Volume=VolumeTetrahedronAndNormalsigns(node); 

x=node(:,1);
y=node(:,2);
z=node(:,3);

q_point(:,1) = x(1) * (1-xw(:,1)-xw(:,2)-xw(:,3)-xw(:,4)) + x(2) * xw(:,1) + x(3) * xw(:,2) + x(4) * xw(:,3);
q_point(:,2) = y(1) * (1-xw(:,1)-xw(:,2)-xw(:,3)-xw(:,4)) + y(2) * xw(:,1) + y(3) * xw(:,2) + y(4) * xw(:,3);
q_point(:,3) = z(1) * (1-xw(:,1)-xw(:,2)-xw(:,3)-xw(:,4)) + z(2) * xw(:,1) + z(3) * xw(:,2) + z(4) * xw(:,3);



number_of_qp=length(q_point(:,1));
[RT0_XX,RT0_XY,RT0_XZ,RT0_YY,RT0_YZ,RT0_ZZ,RT0_div,...
 RT0_XGradX,RT0_XGradY,RT0_XGradZ,RT0_YGradX,RT0_YGradY,RT0_YGradZ,RT0_ZGradX,RT0_ZGradY,RT0_ZGradZ,...
 GradXGradX,GradXGradY,GradXGradZ,GradYGradX,GradYGradY,GradYGradZ,GradZGradX,GradZGradY,GradZGradZ,weights] = phiRT3D(qrule,node);

Valss11=zeros(number_of_qp,1);
Valss12=zeros(number_of_qp,1);
Valss13=zeros(number_of_qp,1);
Valss11=zeros(number_of_qp,1);

 F=zeros(3*face_per_elem,1);

C1=beta*beta;
C2=(2*alpha^2+(alpha+beta)^2);
C3=(alpha^2+2*alpha*(alpha+beta));
C4 = - Cconst*(alpha+beta);
C5 = - 0.5 * Cconst * beta;
C6 = - Cconst * alpha;
C7= Volume * Cconst;
C8=(Cconst * C3 - 2 * Casym);
C9=(Cconst * C1 +2 * Casym);
                     
            % qui se hai f1, f2 devi fare tutti alla fine....
            %for qp=1:number_of_qp
            Valss11= C9 .* (RT0_YY+RT0_ZZ) + Cconst * C2 * RT0_XX + Ceq * RT0_div  ;
                 
            Valss12= C8 * RT0_XY;

            Valss13= C8* RT0_XZ;
                 
                 
            Valss22= C9 * (RT0_XX  + RT0_ZZ) + Cconst * C2 * RT0_YY  + Ceq * RT0_div ;    
                 
            Valss23=C8 * RT0_YZ;  
                 
            Valss33= C9 * RT0_XX + RT0_YY + Cconst * C2 * RT0_ZZ + Ceq * RT0_div ;                 

                 
            Valsu11= (C4 * RT0_XGradX  + C5 *  ( RT0_YGradY + RT0_ZGradZ ))';
                 
            Valsu12= (C6 * RT0_XGradY + C5 * RT0_YGradX)';
                   
            Valsu13= (C6 * RT0_XGradZ + C5 * RT0_ZGradX)';

            Valsu21= (C6 * RT0_YGradX + C5 * RT0_XGradY)';
                 
            Valsu22= (C4 * RT0_YGradY + C5 * ( RT0_XGradX + RT0_ZGradZ ))';
                   
            Valsu23= (C6 * RT0_YGradZ + C5 * RT0_ZGradY)';

            Valsu31= (C6 * RT0_ZGradX + C5 * RT0_XGradZ )';
                 
            Valsu32= (C6 * RT0_ZGradY + C5 * RT0_YGradZ)';
                   
            Valsu33= (C4 * RT0_ZGradZ + C5 * ( RT0_XGradX + RT0_YGradY  ))';
                  
 
            Valuu11= C7 * (    GradXGradX + 0.5 * (GradYGradY+GradZGradZ))' ;  
                 
            Valuu12= C7 * (0.5 * GradXGradY)';

            Valuu13= C7 * (0.5 * GradXGradZ)';
                 
            Valuu22= C7 * (0.5 * GradXGradX +     GradYGradY + 0.5 * GradZGradZ)';  
                 
            Valuu23= C7 * (0.5 * GradYGradZ)';

            Valuu33= C7 * (0.5 * GradXGradX + 0.5 * GradYGradY +        GradZGradZ)';
            
            Valss11=[Valss11(1:4);Valss11(5:8);Valss11(9:12);Valss11(13:16)];
            Valss12=[Valss12(1:4);Valss12(5:8);Valss12(9:12);Valss12(13:16)];
            Valss13=[Valss13(1:4);Valss13(5:8);Valss13(9:12);Valss13(13:16)];
            Valss22=[Valss22(1:4);Valss22(5:8);Valss22(9:12);Valss22(13:16)];
            Valss23=[Valss23(1:4);Valss23(5:8);Valss23(9:12);Valss23(13:16)];
            Valss33=[Valss33(1:4);Valss33(5:8);Valss33(9:12);Valss33(13:16)];

            Valsu11=[Valsu11(1:4);Valsu11(5:8);Valsu11(9:12);Valsu11(13:16)];
            Valsu12=[Valsu12(1:4);Valsu12(5:8);Valsu12(9:12);Valsu12(13:16)];
            Valsu13=[Valsu13(1:4);Valsu13(5:8);Valsu13(9:12);Valsu13(13:16)];
            Valsu21=[Valsu21(1:4);Valsu21(5:8);Valsu21(9:12);Valsu21(13:16)];
            Valsu22=[Valsu22(1:4);Valsu22(5:8);Valsu22(9:12);Valsu22(13:16)];
            Valsu23=[Valsu23(1:4);Valsu23(5:8);Valsu23(9:12);Valsu23(13:16)];
            Valsu31=[Valsu31(1:4);Valsu31(5:8);Valsu31(9:12);Valsu31(13:16)];
            Valsu32=[Valsu32(1:4);Valsu32(5:8);Valsu32(9:12);Valsu32(13:16)];                        
            Valsu33=[Valsu33(1:4);Valsu33(5:8);Valsu33(9:12);Valsu33(13:16)];
            
            
            Valuu11=[Valuu11(1:4);Valuu11(5:8);Valuu11(9:12);Valuu11(13:16)];
            Valuu12=[Valuu12(1:4);Valuu12(5:8);Valuu12(9:12);Valuu12(13:16)];
            Valuu13=[Valuu13(1:4);Valuu13(5:8);Valuu13(9:12);Valuu13(13:16)];
            Valuu22=[Valuu22(1:4);Valuu22(5:8);Valuu22(9:12);Valuu22(13:16)];
            Valuu23=[Valuu23(1:4);Valuu23(5:8);Valuu23(9:12);Valuu23(13:16)];
            Valuu33=[Valuu33(1:4);Valuu33(5:8);Valuu33(9:12);Valuu33(13:16)];

            MS=[Valss11  Valss12 Valss13;
                Valss12' Valss22 Valss23;
                Valss13' Valss23' Valss33;]*Volume;
            
           
            MSU=[Valsu11 Valsu12 Valsu13;
                 Valsu21 Valsu22 Valsu23;
                 Valsu31 Valsu32 Valsu33;]*Volume;
            
            
            MU=[Valuu11  Valuu12 Valuu13;
                Valuu12' Valuu22 Valuu23;
                Valuu13' Valuu23' Valuu33;];
                

   vec=1:face_per_elem;
   
                F(vec,1)= F(vec,1)+ sum( weights .* (-1.0)*Ceq * RT0_div(vec) .* fx(q_point(:,1),q_point(:,2),q_point(:,3)) ,1)';
                F(vec+face_per_elem,1)= F(vec+face_per_elem,1)+ sum( weights .* (-1.0)*Ceq * RT0_div(vec) .* fy(q_point(:,1),q_point(:,2),q_point(:,3)) ,1)';
                F(vec+2*face_per_elem,1)= F(vec+2*face_per_elem,1)+ sum( weights .* (-1.0)*Ceq * RT0_div(vec) .* fz(q_point(:,1),q_point(:,2),q_point(:,3)) ,1)';
  
                
    F= Volume * F;              


          
          
          
end
