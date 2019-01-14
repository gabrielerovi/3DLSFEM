function [M]=assembling_Contact2(C,RT0_mass,RT0_div,P1_Grad,RT0_Grad,RT0_divergence,node,face_per_elem,Ceq)

%global M;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MS = \int_k Asigma : Atau = \int_k Asigma : Atau = \sum_{i=1}^4 \sum_{j=1}^4 \int_k Aphi_i : Aphi_j =                    %%%
%%%%%    = \sum_{i=1}^4 \sum_{j=1}^4 \int_k Aphi_i : Aphi_j = \sum_{i=1}^4 \sum_j=1}^4 |detJ|/(detJ)^2 \int_K AJPhi_i : AJPhi_j%%%
%%%%%    = \sum_{i=1}^4 \sum_j=1}^4 |detJ|/(detJ)^2 \int_K (Beta JPhi_i  +alpha tr(JPhi_i) : (Beta JPhi_j +alpha tr(JPhi_j)    %%%
%%%%%    = \sum_{i=1}^4 \sum_j=1}^4 |detJ|/(detJ)^2 \int_K (Beta JPhi_i  +alpha tr(JPhi_i) : (Beta JPhi_j +alpha tr(JPhi_j)    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MS_{i,j}= |detJ|/(detJ)^2 \int_K (Beta JPhi_i  +alpha tr(JPhi_i) : (Beta JPhi_j +alpha tr(JPhi_j)                      %%%%%  
%%%%% MS_{i,j}_{m,m}=|detJ|/(detJ)^2 \int_K (Beta (Jmx Phi_{i,x}+JmyPhi_{i,y}+Jmz Phi_{i,z})+alpha (Jxx Phi_{i,x}+Jyy Phi_{i,y} + Phi_{i,z}) : (Beta JPhi_j +alpha tr(JPhi_j)

MS=zeros(face_per_elem*3,face_per_elem*3);
MSU=zeros(face_per_elem*3,face_per_elem*3);
MU=zeros(face_per_elem*3,face_per_elem*3);


  [J,detJ]=J_and_detJ(node);   
  absdetJ=abs(detJ);
  Volume=abs(detJ)/6.0;
  frac1_absdetJ=1.0/abs(detJ);
  frac1_detJ=1.0/(detJ);
  frac1_absdetJ2=1.0/abs(detJ)^2;

  
            RT0=RT0_mass*kron(J,J)'*frac1_absdetJ;
            RT0DIV=RT0_div'*frac1_absdetJ;
              
%             qrule=4;
%             [q_point,weights,Volume]=quadrature_points_3D(qrule,node);
%             q_point=q_point';
%             qrule_npoints=length(q_point(1,:));
%             RT0_basisX=repmat(q_point(1,:),face_per_elem,1)-repmat(node(:,1),1,qrule_npoints);
%                 RT0_basisY=repmat(q_point(2,:),face_per_elem,1)-repmat(node(:,2),1,qrule_npoints);
%                 RT0_basisZ=repmat(q_point(3,:),face_per_elem,1)-repmat(node(:,3),1,qrule_npoints);
% 
%               sign_normal=repmat(RT0_divergence,1,qrule_npoints);
%                 RT0_basisX=RT0_basisX.*sign_normal;
%                 RT0_basisY=RT0_basisY.*sign_normal;
%                 RT0_basisZ=RT0_basisZ.*sign_normal;
                
            %C(2)=  C(2)*9.548611111111172e-01  
              
            Valss11= (C(9) .* (RT0(:,5)+RT0(:,9)) + C(2) * RT0(:,1) + Ceq * RT0DIV );%*Volume ;
                             
%             Valss12= (C(8) * RT0(:,2))*Volume;
            Valss12= (C(10) * RT0(:,2)+C(11)*RT0(:,4));%*Volume;

%             Valss13= (C(8)* RT0(:,3))*Volume;
            Valss13= (C(10)* RT0(:,3)+C(11)*RT0(:,7));%*Volume;
                 
            Valss22= (C(9) * (RT0(:,1)  + RT0(:,9)) + C(2) * RT0(:,5)  + Ceq * RT0DIV );%*Volume;    
                 
%             Valss23= (C(8) * RT0(:,6))*Volume;  
            Valss23= (C(10) * RT0(:,6)+C(11)*RT0(:,8));%*Volume;  
            
            Valss33= (C(9) * (RT0(:,1) + RT0(:,5)) + C(2) * RT0(:,9) + Ceq * RT0DIV );%*Volume;                 

            Valss11=reshape(Valss11,[face_per_elem,face_per_elem])';%[Valss11(1:4),Valss11(5:8),Valss11(9:12),Valss11(13:16)];
            Valss12=reshape(Valss12,[face_per_elem,face_per_elem])';
            Valss13=reshape(Valss13,[face_per_elem,face_per_elem])';
            Valss22=reshape(Valss22,[face_per_elem,face_per_elem])';
            Valss23=reshape(Valss23,[face_per_elem,face_per_elem])';
            Valss33=reshape(Valss33,[face_per_elem,face_per_elem])';            
            
            
            MS(1:4,1:4)=Valss11;
            MS(1:4,5:8)=Valss12;
            MS(1:4,9:12)=Valss13;
            MS(5:8,1:4)=Valss12';
            MS(5:8,5:8)=Valss22;
            MS(5:8,9:12)=Valss23;
            MS(9:12,1:4)=Valss13';
            MS(9:12,5:8)=Valss23';
            MS(9:12,9:12)=Valss33;            
            
            
            
      
            
            
            
            
            
            
            
               
            
             JT=frac1_detJ* [J(3,3)*J(2,2)-J(3,2)*J(2,3)        -(J(3,3)*J(2,1)-J(3,1)*J(2,3))  J(3,2)*J(2,1)-J(3,1)*J(2,2);
                             -(J(3,3)*J(1,2)-J(3,2)*J(1,3))  J(3,3)*J(1,1)-J(3,1)*J(1,3)    -(J(3,2)*J(1,1)-J(3,1)*J(1,2));
                             J(2,3)*J(1,2)-J(2,2)*J(1,3)    -(J(2,3)*J(1,1)-J(2,1)*J(1,3))     J(2,2)*J(1,1)-J(2,1)*J(1,2)];        
            %RT0GRAD=RT0_Grad*kron(J,JT)'*frac1_detJ;
            RT0GRAD=RT0_Grad*kron(J,JT)'*sign(detJ);
            
           
            Valsu11= ((C(4) * RT0GRAD(:,1)  + C(5) *  ( RT0GRAD(:,5) + RT0GRAD(:,9) ))');
                 
            Valsu12= ((C(6) * RT0GRAD(:,2) + C(5) * RT0GRAD(:,4))');
                   
            Valsu13= ((C(6) * RT0GRAD(:,3) + C(5) * RT0GRAD(:,7))');

            Valsu21= ((C(6) * RT0GRAD(:,4) + C(5) * RT0GRAD(:,2))');
                 
            Valsu22= ((C(4) * RT0GRAD(:,5) + C(5) * ( RT0GRAD(:,1) + RT0GRAD(:,9) ))');
                   
            Valsu23= ((C(6) * RT0GRAD(:,6) + C(5) * RT0GRAD(:,8))');

            Valsu31= ((C(6) * RT0GRAD(:,7) + C(5) * RT0GRAD(:,3) )');
                 
            Valsu32= ((C(6) * RT0GRAD(:,8) + C(5) * RT0GRAD(:,6))');
                   
            Valsu33= ((C(4) * RT0GRAD(:,9) + C(5) * ( RT0GRAD(:,1) + RT0GRAD(:,5)  ))' );
                  
            Valsu11=reshape(Valsu11,[face_per_elem,face_per_elem])';%[Valsu11(1:4);Valsu11(5:8);Valsu11(9:12);Valsu11(13:16)];
            Valsu12=reshape(Valsu12,[face_per_elem,face_per_elem])';
            Valsu13=reshape(Valsu13,[face_per_elem,face_per_elem])';
            Valsu21=reshape(Valsu21,[face_per_elem,face_per_elem])';
            Valsu22=reshape(Valsu22,[face_per_elem,face_per_elem])';
            Valsu23=reshape(Valsu23,[face_per_elem,face_per_elem])';
            Valsu31=reshape(Valsu31,[face_per_elem,face_per_elem])';
            Valsu32=reshape(Valsu32,[face_per_elem,face_per_elem])';                       
            Valsu33=reshape(Valsu33,[face_per_elem,face_per_elem])';
            
            MSU(1:4,1:4)=Valsu11;
            MSU(1:4,5:8)=Valsu12;
            MSU(1:4,9:12)=Valsu13;
            MSU(5:8,1:4)=Valsu21;
            MSU(5:8,5:8)=Valsu22;
            MSU(5:8,9:12)=Valsu23;
            MSU(9:12,1:4)=Valsu31;
            MSU(9:12,5:8)=Valsu32;
            MSU(9:12,9:12)=Valsu33;
            
            
            
            
            
            
            
            
            GRADGRAD=P1_Grad*kron(JT,JT)';
            C(7)= absdetJ * C(7);
            
            Valuu11= C(7) * (    GRADGRAD(:,1) + 0.5 * (GRADGRAD(:,5)+GRADGRAD(:,9)))' ;  
                 
            %Valuu12= C(7) * (0.5 * GRADGRAD(:,2))';
            Valuu12= C(7) * (0.5 * GRADGRAD(:,4))';

            %Valuu13= C(7) * (0.5 * GRADGRAD(:,3))';
            Valuu13= C(7) * (0.5 * GRADGRAD(:,7))';
                            
            Valuu22= C(7) * (0.5 * GRADGRAD(:,1) +     GRADGRAD(:,5) + 0.5 * GRADGRAD(:,9))';  
                 
            %Valuu23= C(7) * (0.5 * GRADGRAD(:,6))';
            Valuu23= C(7) * (0.5 * GRADGRAD(:,8))';

            Valuu33= C(7) * (0.5 * GRADGRAD(:,1) + 0.5 * GRADGRAD(:,5) +        GRADGRAD(:,9))';
    
  
  
  
                   
 
  
  
        


            
            
            
            Valuu11=reshape(Valuu11,[face_per_elem,face_per_elem])';%[Valuu11(1:4);Valuu11(5:8);Valuu11(9:12);Valuu11(13:16)];
            Valuu12=reshape(Valuu12,[face_per_elem,face_per_elem])';
            Valuu13=reshape(Valuu13,[face_per_elem,face_per_elem])';
            Valuu22=reshape(Valuu22,[face_per_elem,face_per_elem])';
            Valuu23=reshape(Valuu23,[face_per_elem,face_per_elem])';
            Valuu33=reshape(Valuu33,[face_per_elem,face_per_elem])';


            

            
            
            MU(1:4,1:4)=Valuu11;
            MU(1:4,5:8)=Valuu12;
            MU(1:4,9:12)=Valuu13;
            MU(5:8,1:4)=Valuu12';
            MU(5:8,5:8)=Valuu22;
            MU(5:8,9:12)=Valuu23;
            MU(9:12,1:4)=Valuu13';
            MU(9:12,5:8)=Valuu23';
            MU(9:12,9:12)=Valuu33;
            
            
              M=[MS,   MSU;
                 MSU', MU];  
             
              M=M(:);

          
          
          
end
