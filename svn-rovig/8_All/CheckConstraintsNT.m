function [Constraint] = CheckConstraintsNT(mesh,maps,parameters)


gap=parameters.gap;
L=size(mesh);
L=L(1);
grid=mesh{L};
N=grid.N;
NE=grid.NE;
node=grid.node;
edge=grid.edge;
E_bc=grid.E_bc;
normal_node=mesh{L}.normal_node;
normal_edge=mesh{L}.normal_edge;
normal_node_contact=mesh{L}.normal_node_contact;
normal_edge_contact=mesh{L}.normal_edge_contact;
E_contact=grid.E_contact;
N_contact=grid.N_contact;
edge=grid.edge;

CheckConstraintE1=sparse(length(E_contact), 2 * NE + 2 * N);
CheckConstraintN1=sparse(length(N_contact), 2 * NE + 2 * N);
CheckConstraintE2=sparse(length(E_contact), 2 * NE + 2 * N);
CheckConstraintN2=sparse(length(N_contact), 2 * NE + 2 * N);

contE=0;

RhsN1=[];
RhsN2=[];
RhsE1=[];
RhsE2=[];
for ee=E_contact        
        % in the midpoint of each edge, we have:
        % u(midpoint) n_obst(midpoint) <= g
        % we can compute the value in the midpoint as an average between the values of the two adjacent nodes:
        % un(midpoint)=u(midpoint) n_obst(midpoint)
        % since we do not have the value of u in the midpoint, we have to
        % reconstruct it. We use normal and tangent component to do so.
        % Unfortunately the local normals and tangent differ. This mean
        % that, before averaging, we need to use the canonical basis
        %  u_midpoint = 0.5*(H_a'u_a + H_b' u_b)
        % u_n(midpoint)=  n_obst(midpoint) * u_midpoint
        %              =  n_obst(midpoint) * 0.5 * (H_a'u_nt_a + H_b' u_nt_b)
        %              =  (0.5 * n_obst(midpoint) * H_a') u_nt_a + 
        %                 (0.5 * n_obst(midpoint) * H_b') u_nt_b 
        
        contE=contE+1;
        vertices=edge(ee,:);
        
        normal_node_a=normal_node_contact{vertices(1)};
        normal_node_b=normal_node_contact{vertices(2)};
        n_obst=normal_edge_contact{ee}; 
        
        [xn,H_a]=HouseHolderTransformation([0;0],normal_node_a);
        [xn,H_b]=HouseHolderTransformation([0;0],normal_node_b);

        tmp_a=0.5 * n_obst'*H_a';
        tmp_b=0.5 * n_obst'*H_b'; 
        
        CheckConstraintE1(contE, [ 2 * NE + vertices(1),  2 * NE + N + vertices(1)])=[tmp_a(1),tmp_a(2)];
        CheckConstraintE1(contE, [ 2 * NE + vertices(2),  2 * NE + N + vertices(2)])=[tmp_b(1),tmp_b(2)];
        
        RhsE1(contE,1)=gap_function( 0.5 * (node(vertices(1),1) + node(vertices(2),1)), 0.5 * (node(vertices(1),2) + node(vertices(2),2)) );


CheckConstraintE2(contE, ee)=1;
RhsE2(contE,1)=0; 
end

Patch_Edge=maps.Patch_Edge{L};

contN=0;
for nn=N_contact        
    
    
%         % on each vertex, we have:
%         % n_obst(node)' sigma(node) n_node(node) <= 0
%         % we can compute the value in the node as an average between the values of the two adjacent edges:
%         
%         % since we do not have the value of sigma_n in the vertex, we have to
%         % reconstruct it. We use normal and tangent component to do so.
%         % Unfortunately the local normals and tangent differ. This means
%         % that, before averaging, we need to use the canonical basis
%         % sigma_n_(vertex) = 0.5* [sigma_a  (n_obst(vertex) n_a) + sigma_b (n_obst(vertex) n_b)]

        
        contN=contN+1;
        edges=[];
        for ee_loc=Patch_Edge{nn}
            if(length(intersect(nn,edge(ee_loc,:)))>0)
             edges=[edges;ee_loc];   
            end
            edges=intersect(find(E_bc>0),edges);
        end
        
        phidotn1=phi_dot_n(grid,edges(1));
        phidotn2=phi_dot_n(grid,edges(2));
        
        n_obst=normal_node_contact{nn};
        
        % HERE ADD OR NOT phidotn depending on how you perform the
        % HOUSEHOLDER transformation H
        CheckConstraintN1(contN,edges(1))=0.5  * n_obst' * normal_edge{edges(1)};
        CheckConstraintN1(contN,edges(2))=0.5  * n_obst' * normal_edge{edges(2)};
        RhsN1(contN,1)=0;


vertices=node(nn,[1,2]);
CheckConstraintN2(contN,2*NE+nn)=1;
RhsN2(contN,1)=gap_function(vertices(1),vertices(2) );
end


  Constraint.CheckConstraintE1=CheckConstraintE1;
  Constraint.RhsE1=RhsE1;
  Constraint.CheckConstraintN1=CheckConstraintN1;
  Constraint.RhsN1=RhsN1;
  Constraint.CheckConstraintE2=CheckConstraintE2;
  Constraint.RhsE2=RhsE2;
  Constraint.CheckConstraintN2=CheckConstraintN2;
  Constraint.RhsN2=RhsN2;
  
  
  Constraint.WorkingSetE=E_contact;
  Constraint.WorkingSetN=[];
  
end