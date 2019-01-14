function [V,sign_normal,Area]=VolumeTetrahedronAndNormalsigns(nodes_coord_of_tetrahedron) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Given a tetrahedron of nodes n1,n2,n3,n4, consider one of its face:                            %%%%
%%%% n2,n3,n4: normal_234 = cross(p3-p2,p4-p2), sgn= sign((p2-p1) * normal_234)                     %%%%
%%%% n1,n3,n4: normal_134 = cross(p3-p1,p4-p1), sgn= sign((p1-p2) * normal_134)                     %%%%
%%%% n1,n2,n4: normal_124 = cross(p2-p1,p4-p1), sgn= sign((p1-p3) * normal_124)                     %%%%
%%%% n1,n2,n3: normal_123 = cross(p2-p1,p3-p1), sgn= sign((p1-p4) * normal_123)                     %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   s21=nodes_coord_of_tetrahedron(2,:)-nodes_coord_of_tetrahedron(1,:);
   s31=nodes_coord_of_tetrahedron(3,:)-nodes_coord_of_tetrahedron(1,:);
   s41=nodes_coord_of_tetrahedron(4,:)-nodes_coord_of_tetrahedron(1,:);
   s32=nodes_coord_of_tetrahedron(3,:)-nodes_coord_of_tetrahedron(2,:);
   s42=nodes_coord_of_tetrahedron(4,:)-nodes_coord_of_tetrahedron(2,:); 
   s43=nodes_coord_of_tetrahedron(4,:)-nodes_coord_of_tetrahedron(3,:);

   normal_234=cross(s32,s42)';
   normal_134 = cross(s31,s41)';
   normal_124 = cross(s21,s41)';
   normal_123 = cross(s21,s31)'; 


   V= abs(s21*cross(s31,s41 )')/6.0;
   sign_normal=sign([s21*normal_234; (-s21)*normal_134; (-s31)*normal_124; (-s41)*normal_123;]); 

   
   for ii=1:4
   face=nodes_coord_of_tetrahedron;
   face(ii,:)=[];
   Area(ii,1)=AreaTriangle(face) ;
   end






end
