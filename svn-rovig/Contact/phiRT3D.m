function [RT0_XX,RT0_XY,RT0_XZ,RT0_YX,RT0_YY,RT0_YZ,RT0_ZX,RT0_ZY,RT0_ZZ,RT0_div,RT0_divergence,signRT0_divergence,...
          RT0_XGradX,RT0_XGradY,RT0_XGradZ,RT0_YGradX,RT0_YGradY,RT0_YGradZ,RT0_ZGradX,RT0_ZGradY,RT0_ZGradZ,...
          GradXGradX,GradXGradY,GradXGradZ,GradYGradX,GradYGradY,GradYGradZ,GradZGradX,GradZGradY,GradZGradZ,q_point,weights] = phiRT3D(qrule,nodes_coord_of_tetrahedron)
% q_point is a matrix m x n
% m = qrule_npoints
% n = dimension of the problem

% node contains the node of the triangle
    [q_point,weights,Volume]=quadrature_points_3D(qrule,nodes_coord_of_tetrahedron);

    q_point=q_point';
    dim=3;
    face_per_elem=dim+1;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Given a tetrahedron of nodes n1,n2,n3,n4, consider one of its face:                            %%%%
%%%% n2,n3,n4: normal_234 = cross(p3-p2,p4-p2), sgn= sign((p2-p1) * normal_234)                     %%%%
%%%% n1,n3,n4: normal_134 = cross(p3-p1,p4-p1), sgn= sign((p1-p2) * normal_134)                     %%%%
%%%% n1,n2,n4: normal_124 = cross(p2-p1,p4-p1), sgn= sign((p1-p3) * normal_124)                     %%%%
%%%% n1,n2,n3: normal_124 = cross(p2-p1,p3-p1), sgn= sign((p1-p3) * normal_123)                     %%%%
%%%% V = abs((p2-p1)'[(p3-p1) x (p4-p1)])/6 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    [Volume,sign_normal,FaceArea]=VolumeTetrahedronAndNormalsigns(nodes_coord_of_tetrahedron) ;
    
    %%%% FENICS frac_1_dimV
    %frac_1_dimV=1.0/(2*dim*Volume);  
    %%%% MY frac_1_dimV
    frac_1_dimV=1.0/(dim*Volume);  
    
%     frac_1_H=FaceArea/(dim*Volume); 
%     a=nodes_coord_of_tetrahedron(1,:);
%     b=nodes_coord_of_tetrahedron(2,:); 
%     c=nodes_coord_of_tetrahedron(3,:); 
%     d=nodes_coord_of_tetrahedron(3,:); 

    qrule_npoints=length(q_point(1,:));


    RT0_basis=zeros(face_per_elem,qrule_npoints,dim);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% sss = # basis function                                         %%%%
        %%%% qp  = quadrature point                                         %%%%
        %%%% ddd = component                                                %%%%
        %%%% Compute phi_i=(x-p_i)/(3 V), with p_i opposite node to face i  %%%% 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                RT0_basisX=repmat(q_point(1,:),face_per_elem,1)-repmat(nodes_coord_of_tetrahedron(:,1),1,qrule_npoints);
                RT0_basisY=repmat(q_point(2,:),face_per_elem,1)-repmat(nodes_coord_of_tetrahedron(:,2),1,qrule_npoints);
                RT0_basisZ=repmat(q_point(3,:),face_per_elem,1)-repmat(nodes_coord_of_tetrahedron(:,3),1,qrule_npoints);
 
                
                
                % ATTENZIONE QUI CAMBIO SEGNO PERCHE' MI SA CHE LA
                % CONVEZIONE DI FENICS E' DIFFERENTE
            %    sign_normal=-sign_normal;
                
                
                RT0_divergence= sign_normal*frac_1_dimV;
                %RT0_divergence= sign_normal.*frac_1_H;
                
                
                

                
                sign_normal=repmat(RT0_divergence,1,qrule_npoints);
                
                
                
                
                
                
                
                
                
                
                

                
                RT0_basisX=RT0_basisX.*sign_normal;
                RT0_basisY=RT0_basisY.*sign_normal;
                RT0_basisZ=RT0_basisZ.*sign_normal;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%                             X- COMPONENTS                                %%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
                RT0_X(1)=sum( RT0_basisX(1,:).*weights');
                RT0_X(2)=sum( RT0_basisX(2,:).*weights');
                RT0_X(3)=sum( RT0_basisX(3,:).*weights');
                RT0_X(4)=sum( RT0_basisX(4,:).*weights');
                
     
              
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%                             Y- COMPONENTS                                %%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
                RT0_Y(1)=sum( RT0_basisY(1,:).*weights');
                RT0_Y(2)=sum( RT0_basisY(2,:).*weights');
                RT0_Y(3)=sum( RT0_basisY(3,:).*weights');
                RT0_Y(4)=sum( RT0_basisY(4,:).*weights');
                
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%                             Z- COMPONENTS                                %%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
                RT0_Z(1)=sum( RT0_basisZ(1,:).*weights');
                RT0_Z(2)=sum( RT0_basisZ(2,:).*weights');
                RT0_Z(3)=sum( RT0_basisZ(3,:).*weights');
                RT0_Z(4)=sum( RT0_basisZ(4,:).*weights');
                               
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%                             XX- COMPONENTS                                %%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
                RT0_XX(1)=sum( (RT0_basisX(1,:).*RT0_basisX(1,:)).*weights');
                RT0_XX(2)=sum( (RT0_basisX(1,:).*RT0_basisX(2,:)).*weights');
                RT0_XX(3)=sum( (RT0_basisX(1,:).*RT0_basisX(3,:)).*weights');
                RT0_XX(4)=sum( (RT0_basisX(1,:).*RT0_basisX(4,:)).*weights');
                
                RT0_XX(5)=sum( (RT0_basisX(2,:).*RT0_basisX(1,:)).*weights');
                RT0_XX(6)=sum( (RT0_basisX(2,:).*RT0_basisX(2,:)).*weights');
                RT0_XX(7)=sum( (RT0_basisX(2,:).*RT0_basisX(3,:)).*weights');
                RT0_XX(8)=sum( (RT0_basisX(2,:).*RT0_basisX(4,:)).*weights');                
                
                RT0_XX(9)=sum( (RT0_basisX(3,:).*RT0_basisX(1,:)).*weights');
                RT0_XX(10)=sum( (RT0_basisX(3,:).*RT0_basisX(2,:)).*weights');
                RT0_XX(11)=sum( (RT0_basisX(3,:).*RT0_basisX(3,:)).*weights');
                RT0_XX(12)=sum( (RT0_basisX(3,:).*RT0_basisX(4,:)).*weights');     
                
                RT0_XX(13)=sum( (RT0_basisX(4,:).*RT0_basisX(1,:)).*weights');
                RT0_XX(14)=sum( (RT0_basisX(4,:).*RT0_basisX(2,:)).*weights');
                RT0_XX(15)=sum( (RT0_basisX(4,:).*RT0_basisX(3,:)).*weights');
                RT0_XX(16)=sum( (RT0_basisX(4,:).*RT0_basisX(4,:)).*weights');     
                
                
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%                             XY- COMPONENTS                                %%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
                RT0_XY(1)=sum( (RT0_basisX(1,:).*RT0_basisY(1,:)).*weights');
                RT0_XY(2)=sum( (RT0_basisX(1,:).*RT0_basisY(2,:)).*weights');
                RT0_XY(3)=sum( (RT0_basisX(1,:).*RT0_basisY(3,:)).*weights');
                RT0_XY(4)=sum( (RT0_basisX(1,:).*RT0_basisY(4,:)).*weights');
                
                RT0_XY(5)=sum( (RT0_basisX(2,:).*RT0_basisY(1,:)).*weights');
                RT0_XY(6)=sum( (RT0_basisX(2,:).*RT0_basisY(2,:)).*weights');
                RT0_XY(7)=sum( (RT0_basisX(2,:).*RT0_basisY(3,:)).*weights');
                RT0_XY(8)=sum( (RT0_basisX(2,:).*RT0_basisY(4,:)).*weights');   
                
                RT0_XY(9)=sum( (RT0_basisX(3,:).*RT0_basisY(1,:)).*weights');
                RT0_XY(10)=sum( (RT0_basisX(3,:).*RT0_basisY(2,:)).*weights');
                RT0_XY(11)=sum( (RT0_basisX(3,:).*RT0_basisY(3,:)).*weights');
                RT0_XY(12)=sum( (RT0_basisX(3,:).*RT0_basisY(4,:)).*weights');   
                
                RT0_XY(13)=sum( (RT0_basisX(4,:).*RT0_basisY(1,:)).*weights');
                RT0_XY(14)=sum( (RT0_basisX(4,:).*RT0_basisY(2,:)).*weights');
                RT0_XY(15)=sum( (RT0_basisX(4,:).*RT0_basisY(3,:)).*weights');            
                RT0_XY(16)=sum( (RT0_basisX(4,:).*RT0_basisY(4,:)).*weights'); 
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%                             XZ- COMPONENTS                                %%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
                RT0_XZ(1)=sum( (RT0_basisX(1,:).*RT0_basisZ(1,:)).*weights');
                RT0_XZ(2)=sum( (RT0_basisX(1,:).*RT0_basisZ(2,:)).*weights');
                RT0_XZ(3)=sum( (RT0_basisX(1,:).*RT0_basisZ(3,:)).*weights');
                RT0_XZ(4)=sum( (RT0_basisX(1,:).*RT0_basisZ(4,:)).*weights');
                
                RT0_XZ(5)=sum( (RT0_basisX(2,:).*RT0_basisZ(1,:)).*weights');
                RT0_XZ(6)=sum( (RT0_basisX(2,:).*RT0_basisZ(2,:)).*weights');
                RT0_XZ(7)=sum( (RT0_basisX(2,:).*RT0_basisZ(3,:)).*weights');
                RT0_XZ(8)=sum( (RT0_basisX(2,:).*RT0_basisZ(4,:)).*weights');  
                
                RT0_XZ(9)=sum( (RT0_basisX(3,:).*RT0_basisZ(1,:)).*weights');
                RT0_XZ(10)=sum( (RT0_basisX(3,:).*RT0_basisZ(2,:)).*weights');            
                RT0_XZ(11)=sum( (RT0_basisX(3,:).*RT0_basisZ(3,:)).*weights');
                RT0_XZ(12)=sum( (RT0_basisX(3,:).*RT0_basisZ(4,:)).*weights'); 
                
                RT0_XZ(13)=sum( (RT0_basisX(4,:).*RT0_basisZ(1,:)).*weights');
                RT0_XZ(14)=sum( (RT0_basisX(4,:).*RT0_basisZ(2,:)).*weights'); 
                RT0_XZ(15)=sum( (RT0_basisX(4,:).*RT0_basisZ(3,:)).*weights'); 
                RT0_XZ(16)=sum( (RT0_basisX(4,:).*RT0_basisZ(4,:)).*weights'); 
                
                
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%                             YX- COMPONENTS                                %%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
                RT0_YX(1)=sum( (RT0_basisY(1,:).*RT0_basisX(1,:)).*weights');
                RT0_YX(2)=sum( (RT0_basisY(1,:).*RT0_basisX(2,:)).*weights');
                RT0_YX(3)=sum( (RT0_basisY(1,:).*RT0_basisX(3,:)).*weights');
                RT0_YX(4)=sum( (RT0_basisY(1,:).*RT0_basisX(4,:)).*weights');
                
                RT0_YX(5)=sum( (RT0_basisY(2,:).*RT0_basisX(1,:)).*weights');
                RT0_YX(6)=sum( (RT0_basisY(2,:).*RT0_basisX(2,:)).*weights');
                RT0_YX(7)=sum( (RT0_basisY(2,:).*RT0_basisX(3,:)).*weights');
                RT0_YX(8)=sum( (RT0_basisY(2,:).*RT0_basisX(4,:)).*weights');  
                
                RT0_YX(9)=sum( (RT0_basisY(3,:).*RT0_basisX(1,:)).*weights');
                RT0_YX(10)=sum( (RT0_basisY(3,:).*RT0_basisX(2,:)).*weights');              
                RT0_YX(11)=sum( (RT0_basisY(3,:).*RT0_basisX(3,:)).*weights');
                RT0_YX(12)=sum( (RT0_basisY(3,:).*RT0_basisX(4,:)).*weights'); 
                
                RT0_YX(13)=sum( (RT0_basisY(4,:).*RT0_basisX(1,:)).*weights');
                RT0_YX(14)=sum( (RT0_basisY(4,:).*RT0_basisX(2,:)).*weights');
                RT0_YX(15)=sum( (RT0_basisY(4,:).*RT0_basisX(3,:)).*weights');
                RT0_YX(16)=sum( (RT0_basisY(4,:).*RT0_basisX(4,:)).*weights');                
                
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%                             YY- COMPONENTS                                %%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
                RT0_YY(1)=sum( (RT0_basisY(1,:).*RT0_basisY(1,:)).*weights');
                RT0_YY(2)=sum( (RT0_basisY(1,:).*RT0_basisY(2,:)).*weights');
                RT0_YY(3)=sum( (RT0_basisY(1,:).*RT0_basisY(3,:)).*weights');
                RT0_YY(4)=sum( (RT0_basisY(1,:).*RT0_basisY(4,:)).*weights');
                
                RT0_YY(5)=sum( (RT0_basisY(2,:).*RT0_basisY(1,:)).*weights');
                RT0_YY(6)=sum( (RT0_basisY(2,:).*RT0_basisY(2,:)).*weights');
                RT0_YY(7)=sum( (RT0_basisY(2,:).*RT0_basisY(3,:)).*weights');
                RT0_YY(8)=sum( (RT0_basisY(2,:).*RT0_basisY(4,:)).*weights');   
                
                
                RT0_YY(9)=sum( (RT0_basisY(3,:).*RT0_basisY(1,:)).*weights');
                RT0_YY(10)=sum( (RT0_basisY(3,:).*RT0_basisY(2,:)).*weights');       
                RT0_YY(11)=sum( (RT0_basisY(3,:).*RT0_basisY(3,:)).*weights');
                RT0_YY(12)=sum( (RT0_basisY(3,:).*RT0_basisY(4,:)).*weights');    
                
                
                RT0_YY(13)=sum( (RT0_basisY(4,:).*RT0_basisY(1,:)).*weights');
                RT0_YY(14)=sum( (RT0_basisY(4,:).*RT0_basisY(2,:)).*weights');
                RT0_YY(15)=sum( (RT0_basisY(4,:).*RT0_basisY(3,:)).*weights');
                RT0_YY(16)=sum( (RT0_basisY(4,:).*RT0_basisY(4,:)).*weights'); 
                

                
                
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%                             YZ- COMPONENTS                                %%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
                RT0_YZ(1)=sum( (RT0_basisY(1,:).*RT0_basisZ(1,:)).*weights');
                RT0_YZ(2)=sum( (RT0_basisY(1,:).*RT0_basisZ(2,:)).*weights');
                RT0_YZ(3)=sum( (RT0_basisY(1,:).*RT0_basisZ(3,:)).*weights');
                RT0_YZ(4)=sum( (RT0_basisY(1,:).*RT0_basisZ(4,:)).*weights');
                
                RT0_YZ(5)=sum( (RT0_basisY(2,:).*RT0_basisZ(1,:)).*weights');          
                RT0_YZ(6)=sum( (RT0_basisY(2,:).*RT0_basisZ(2,:)).*weights');
                RT0_YZ(7)=sum( (RT0_basisY(2,:).*RT0_basisZ(3,:)).*weights');
                RT0_YZ(8)=sum( (RT0_basisY(2,:).*RT0_basisZ(4,:)).*weights'); 
                
                RT0_YZ(9)=sum( (RT0_basisY(3,:).*RT0_basisZ(1,:)).*weights');
                RT0_YZ(10)=sum( (RT0_basisY(3,:).*RT0_basisZ(2,:)).*weights');          
                RT0_YZ(11)=sum( (RT0_basisY(3,:).*RT0_basisZ(3,:)).*weights');
                RT0_YZ(12)=sum( (RT0_basisY(3,:).*RT0_basisZ(4,:)).*weights'); 
                
                RT0_YZ(13)=sum( (RT0_basisY(4,:).*RT0_basisZ(1,:)).*weights');
                RT0_YZ(14)=sum( (RT0_basisY(4,:).*RT0_basisZ(2,:)).*weights');
                RT0_YZ(15)=sum( (RT0_basisY(4,:).*RT0_basisZ(3,:)).*weights');
                RT0_YZ(16)=sum( (RT0_basisY(4,:).*RT0_basisZ(4,:)).*weights'); 

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%                             ZX- COMPONENTS                                %%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
                RT0_ZX(1)=sum( (RT0_basisZ(1,:).*RT0_basisX(1,:)).*weights');
                RT0_ZX(2)=sum( (RT0_basisZ(1,:).*RT0_basisX(2,:)).*weights');
                RT0_ZX(3)=sum( (RT0_basisZ(1,:).*RT0_basisX(3,:)).*weights');
                RT0_ZX(4)=sum( (RT0_basisZ(1,:).*RT0_basisX(4,:)).*weights');
                
                RT0_ZX(5)=sum( (RT0_basisZ(2,:).*RT0_basisX(1,:)).*weights');
                RT0_ZX(6)=sum( (RT0_basisZ(2,:).*RT0_basisX(2,:)).*weights');
                RT0_ZX(7)=sum( (RT0_basisZ(2,:).*RT0_basisX(3,:)).*weights');
                RT0_ZX(8)=sum( (RT0_basisZ(2,:).*RT0_basisX(4,:)).*weights');  
                
                RT0_ZX(9)=sum( (RT0_basisZ(3,:).*RT0_basisX(1,:)).*weights');
                RT0_ZX(10)=sum( (RT0_basisZ(3,:).*RT0_basisX(2,:)).*weights');            
                RT0_ZX(11)=sum( (RT0_basisZ(3,:).*RT0_basisX(3,:)).*weights');
                RT0_ZX(12)=sum( (RT0_basisZ(3,:).*RT0_basisX(4,:)).*weights'); 
                
                RT0_ZX(13)=sum( (RT0_basisZ(4,:).*RT0_basisX(1,:)).*weights');
                RT0_ZX(14)=sum( (RT0_basisZ(4,:).*RT0_basisX(2,:)).*weights');
                RT0_ZX(15)=sum( (RT0_basisZ(4,:).*RT0_basisX(3,:)).*weights');
                RT0_ZX(16)=sum( (RT0_basisZ(4,:).*RT0_basisX(4,:)).*weights'); 

                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%                             ZY- COMPONENTS                                %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
                RT0_ZY(1)=sum( (RT0_basisZ(1,:).*RT0_basisY(1,:)).*weights');
                RT0_ZY(2)=sum( (RT0_basisZ(1,:).*RT0_basisY(2,:)).*weights');
                RT0_ZY(3)=sum( (RT0_basisZ(1,:).*RT0_basisY(3,:)).*weights');
                RT0_ZY(4)=sum( (RT0_basisZ(1,:).*RT0_basisY(4,:)).*weights');
                
                RT0_ZY(5)=sum( (RT0_basisZ(2,:).*RT0_basisY(1,:)).*weights');
                RT0_ZY(6)=sum( (RT0_basisZ(2,:).*RT0_basisY(2,:)).*weights');
                RT0_ZY(7)=sum( (RT0_basisZ(2,:).*RT0_basisY(3,:)).*weights');
                RT0_ZY(8)=sum( (RT0_basisZ(2,:).*RT0_basisY(4,:)).*weights');  
                
                RT0_ZY(9)=sum( (RT0_basisZ(3,:).*RT0_basisY(1,:)).*weights');
                RT0_ZY(10)=sum( (RT0_basisZ(3,:).*RT0_basisY(2,:)).*weights');              
                RT0_ZY(11)=sum( (RT0_basisZ(3,:).*RT0_basisY(3,:)).*weights');
                RT0_ZY(12)=sum( (RT0_basisZ(3,:).*RT0_basisY(4,:)).*weights'); 
                
                RT0_ZY(13)=sum( (RT0_basisZ(4,:).*RT0_basisY(1,:)).*weights');
                RT0_ZY(14)=sum( (RT0_basisZ(4,:).*RT0_basisY(2,:)).*weights');
                RT0_ZY(15)=sum( (RT0_basisZ(4,:).*RT0_basisY(3,:)).*weights'); 
                RT0_ZY(16)=sum( (RT0_basisZ(4,:).*RT0_basisY(4,:)).*weights'); 

                
                
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%                             ZZ- COMPONENTS                                %%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
                RT0_ZZ(1)=sum( (RT0_basisZ(1,:).*RT0_basisZ(1,:)).*weights');
                RT0_ZZ(2)=sum( (RT0_basisZ(1,:).*RT0_basisZ(2,:)).*weights');
                RT0_ZZ(3)=sum( (RT0_basisZ(1,:).*RT0_basisZ(3,:)).*weights');
                RT0_ZZ(4)=sum( (RT0_basisZ(1,:).*RT0_basisZ(4,:)).*weights');
                
                RT0_ZZ(5)=sum( (RT0_basisZ(2,:).*RT0_basisZ(1,:)).*weights');
                RT0_ZZ(6)=sum( (RT0_basisZ(2,:).*RT0_basisZ(2,:)).*weights');
                RT0_ZZ(7)=sum( (RT0_basisZ(2,:).*RT0_basisZ(3,:)).*weights');
                RT0_ZZ(8)=sum( (RT0_basisZ(2,:).*RT0_basisZ(4,:)).*weights');     
                
                RT0_ZZ(9)=sum( (RT0_basisZ(3,:).*RT0_basisZ(1,:)).*weights');
                RT0_ZZ(10)=sum( (RT0_basisZ(3,:).*RT0_basisZ(2,:)).*weights');                  
                RT0_ZZ(11)=sum( (RT0_basisZ(3,:).*RT0_basisZ(3,:)).*weights');
                RT0_ZZ(12)=sum( (RT0_basisZ(3,:).*RT0_basisZ(4,:)).*weights');  
                
                RT0_ZZ(13)=sum( (RT0_basisZ(4,:).*RT0_basisZ(1,:)).*weights');
                RT0_ZZ(14)=sum( (RT0_basisZ(4,:).*RT0_basisZ(2,:)).*weights');
                RT0_ZZ(15)=sum( (RT0_basisZ(4,:).*RT0_basisZ(3,:)).*weights');             
                RT0_ZZ(16)=sum( (RT0_basisZ(4,:).*RT0_basisZ(4,:)).*weights');               
                
                
                
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%                             DIV- COMPONENTS                               %%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              

                RT0_divergence= RT0_divergence * dim;
                signRT0_divergence=sign(RT0_divergence);
                RT0_div(1)=RT0_divergence(1)*RT0_divergence(1);
                RT0_div(2)=RT0_divergence(1)*RT0_divergence(2);
                RT0_div(3)=RT0_divergence(1)*RT0_divergence(3);
                RT0_div(4)=RT0_divergence(1)*RT0_divergence(4);
                
                RT0_div(5)=RT0_div(2);
                RT0_div(6)=RT0_divergence(2)*RT0_divergence(2);
                RT0_div(7)=RT0_divergence(2)*RT0_divergence(3);
                RT0_div(8)=RT0_divergence(2)*RT0_divergence(4);

                RT0_div(9)=RT0_div(3);
                RT0_div(10)=RT0_div(7);
                RT0_div(11)=RT0_divergence(3)*RT0_divergence(3);
                RT0_div(12)=RT0_divergence(3)*RT0_divergence(4);
                
                RT0_div(13)=RT0_div(4);
                RT0_div(14)=RT0_div(8);
                RT0_div(15)=RT0_div(12);
                RT0_div(16)=RT0_divergence(4)*RT0_divergence(4);
                
                RT0_div=RT0_div*Volume;
                
                
                
                
                
                
               
                [P1grad] = P1grad3D(nodes_coord_of_tetrahedron);

                RT0_XGradX=[P1grad(:,1).*RT0_X]; RT0_XGradX=RT0_XGradX(:);
                RT0_XGradY=[P1grad(:,2).*RT0_X]; RT0_XGradY=RT0_XGradY(:);
                RT0_XGradZ=[P1grad(:,3).*RT0_X]; RT0_XGradZ=RT0_XGradZ(:);
                
                RT0_YGradX=[P1grad(:,1).*RT0_Y]; RT0_YGradX=RT0_YGradX(:);
                RT0_YGradY=[P1grad(:,2).*RT0_Y]; RT0_YGradY=RT0_YGradY(:);
                RT0_YGradZ=[P1grad(:,3).*RT0_Y]; RT0_YGradZ=RT0_YGradZ(:);
                
                RT0_ZGradX=[P1grad(:,1).*RT0_Z]; RT0_ZGradX=RT0_ZGradX(:);
                RT0_ZGradY=[P1grad(:,2).*RT0_Z]; RT0_ZGradY=RT0_ZGradY(:);
                RT0_ZGradZ=[P1grad(:,3).*RT0_Z]; RT0_ZGradZ=RT0_ZGradZ(:);                
                

%                 RT0_XGradX=RT0_XGradX*Volume;
%                 RT0_XGradY=RT0_XGradY*Volume;
%                 RT0_XGradZ=RT0_XGradZ*Volume;
                
%                 RT0_YGradX=RT0_YGradX*Volume;
%                 RT0_YGradY=RT0_YGradY*Volume;
%                 RT0_YGradZ=RT0_YGradZ*Volume;
%                 
%                 RT0_ZGradX=RT0_ZGradX*Volume;
%                 RT0_ZGradY=RT0_ZGradY*Volume;
%                 RT0_ZGradZ=RT0_ZGradZ*Volume;
                
                
                P1gradX=repmat(P1grad(:,1)',4,1);
                P1gradY=repmat(P1grad(:,2)',4,1);
                P1gradZ=repmat(P1grad(:,3)',4,1);                

                GradXGradX=[P1grad(:,1)*P1grad(1,1);P1grad(:,1)*P1grad(2,1);P1grad(:,1)*P1grad(3,1);P1grad(:,1)*P1grad(4,1)];
                GradYGradX=[P1grad(:,1)*P1grad(1,2);P1grad(:,1)*P1grad(2,2);P1grad(:,1)*P1grad(3,2);P1grad(:,1)*P1grad(4,2)];
                GradZGradX=[P1grad(:,1)*P1grad(1,3);P1grad(:,1)*P1grad(2,3);P1grad(:,1)*P1grad(3,3);P1grad(:,1)*P1grad(4,3)];

                GradXGradY=[P1grad(:,2)*P1grad(1,1);P1grad(:,2)*P1grad(2,1);P1grad(:,2)*P1grad(3,1);P1grad(:,2)*P1grad(4,1)];
                GradYGradY=[P1grad(:,2)*P1grad(1,2);P1grad(:,2)*P1grad(2,2);P1grad(:,2)*P1grad(3,2);P1grad(:,2)*P1grad(4,2)];
                GradZGradY=[P1grad(:,2)*P1grad(1,3);P1grad(:,2)*P1grad(2,3);P1grad(:,2)*P1grad(3,3);P1grad(:,2)*P1grad(4,3)];

                GradXGradZ=[P1grad(:,3)*P1grad(1,1);P1grad(:,3)*P1grad(2,1);P1grad(:,3)*P1grad(3,1);P1grad(:,3)*P1grad(4,1)];
                GradYGradZ=[P1grad(:,3)*P1grad(1,2);P1grad(:,3)*P1grad(2,2);P1grad(:,3)*P1grad(3,2);P1grad(:,3)*P1grad(4,2)];
                GradZGradZ=[P1grad(:,3)*P1grad(1,3);P1grad(:,3)*P1grad(2,3);P1grad(:,3)*P1grad(3,3);P1grad(:,3)*P1grad(4,3)];

                
                GradXGradX=GradXGradX*Volume;
                GradXGradY=GradXGradY*Volume;
                GradXGradZ=GradXGradZ*Volume;
                
                GradYGradX=GradYGradX*Volume;
                GradYGradY=GradYGradY*Volume;
                GradYGradZ=GradYGradZ*Volume;
                
                GradZGradX=GradZGradX*Volume;
                GradZGradY=GradZGradY*Volume;
                GradZGradZ=GradZGradZ*Volume;                
                
                
                
                
end
        
        