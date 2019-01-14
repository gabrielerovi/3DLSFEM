function [RT0_basis,RT0_divergence] = phiRT3Dcell(q_point,nodes_coord_of_tetrahedron)
% q_point is a matrix m x n
% m = qrule_npoints
% n = dimension of the problem

% node contains the node of the triangle
    
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
   
    [Volume,sign_normal]=VolumeTetrahedronAndNormalsigns(nodes_coord_of_tetrahedron) ;
    frac_1_dimV=1.0/(dim*Volume);    
%     a=nodes_coord_of_tetrahedron(1,:);
%     b=nodes_coord_of_tetrahedron(2,:); 
%     c=nodes_coord_of_tetrahedron(3,:); 
%     d=nodes_coord_of_tetrahedron(3,:); 

    qrule_npoints=length(q_point(:,1));


    RT0_basis=sparse(face_per_elem,qrule_npoints,dim);

    RT0_basis=cell(face_per_elem,1);
    for qp=1:qrule_npoints
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% sss = # basis function                                         %%%%
        %%%% qp  = quadrature point                                         %%%%
        %%%% ddd = component                                                %%%%
        %%%% Compute phi_i=(x-p_i)/(3 V), with p_i opposite node to face i  %%%% 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for ddd=1:dim
            for sss=1:face_per_elem
                RT0_basis{sss}(qp,ddd)=sign_normal(sss) * (q_point(qp,ddd)-nodes_coord_of_tetrahedron(sss,ddd))*frac_1_dimV ;
            end
        end

     end

        for sss=1:face_per_elem
        RT0_divergence(sss)= sign_normal(sss) * dim*frac_1_dimV ;  
        end
        
end
        
        