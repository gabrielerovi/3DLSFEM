function [F]=assembling_ContactRhs(RT0_divergence,node,face_per_elem,Ceq,fx,fy,fz,q_pointref,weights,istheexternalforcenonzero)

%global M;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MS = \int_k Asigma : Atau = \int_k Asigma : Atau = \sum_{i=1}^4 \sum_{j=1}^4 \int_k Aphi_i : Aphi_j =                    %%%
%%%%%    = \sum_{i=1}^4 \sum_{j=1}^4 \int_k Aphi_i : Aphi_j = \sum_{i=1}^4 \sum_j=1}^4 |detJ|/(detJ)^2 \int_K AJPhi_i : AJPhi_j%%%
%%%%%    = \sum_{i=1}^4 \sum_j=1}^4 |detJ|/(detJ)^2 \int_K (Beta JPhi_i  +alpha tr(JPhi_i) : (Beta JPhi_j +alpha tr(JPhi_j)    %%%
%%%%%    = \sum_{i=1}^4 \sum_j=1}^4 |detJ|/(detJ)^2 \int_K (Beta JPhi_i  +alpha tr(JPhi_i) : (Beta JPhi_j +alpha tr(JPhi_j)    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MS_{i,j}= |detJ|/(detJ)^2 \int_K (Beta JPhi_i  +alpha tr(JPhi_i) : (Beta JPhi_j +alpha tr(JPhi_j)                      %%%%%  
%%%%% MS_{i,j}_{m,m}=|detJ|/(detJ)^2 \int_K (Beta (Jmx Phi_{i,x}+JmyPhi_{i,y}+Jmz Phi_{i,z})+alpha (Jxx Phi_{i,x}+Jyy Phi_{i,y} + Phi_{i,z}) : (Beta JPhi_j +alpha tr(JPhi_j)




  [J,detJ]=J_and_detJ(node);   
  absdetJ=abs(detJ);
  Volume=abs(detJ)/6.0;
  frac1_absdetJ=1.0/abs(detJ);
  frac1_detJ=1.0/(detJ);
  frac1_absdetJ2=1.0/abs(detJ)^2;

  
              
              
             
     
  
  
  
                   
          
            
                

 
  F=zeros(face_per_elem*3,1);
  % here it should be localRT0DIV*weights*qpoints * Volume. Instead of defining
  % localRT0DIV=sign/Volume, we maintain sign.
   q_point=q_pointref*J'+repmat(node(1,:),length(q_pointref),1);
   
   %%%%% div tau= alpha dim/(dim V)= alpha /V 
   %%%%% V=detJ/6
   %%%%% int_V div tau f= alpha/V int_V f=(alpha/V) detJ weight*f(qpoint)=
   %%%%% alpha 6  weight*f(qpoint)
  vec=1:face_per_elem;
  % 3= divergence on the reference element
                F(vec,1)                 = F(vec,1)                - sign(detJ)*Ceq * RT0_divergence * weights' * fx(q_point(:,1),q_point(:,2),q_point(:,3));
                F(vec+face_per_elem,1)   = F(vec+face_per_elem,1)  - sign(detJ)*Ceq * RT0_divergence * weights' * fy(q_point(:,1),q_point(:,2),q_point(:,3));
                F(vec+2*face_per_elem,1) = F(vec+2*face_per_elem,1)- sign(detJ)*Ceq * RT0_divergence * weights' * fz(q_point(:,1),q_point(:,2),q_point(:,3)) ;


          
          
          
end
