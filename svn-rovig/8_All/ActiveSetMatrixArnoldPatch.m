function [A,b,ST,STT,CT,B,BT,C] = ActiveSetMatrixArnoldPatch(mesh,parameters,A,b,N_contact,E_contact,NE_loc,N_loc)
% we solve for min H, with H=0.5 x' A x - x' f - lambda (B x -c)
% structure of the problem
% |A -B'| |x     |= |f|
% |B  0 | |lambda|  |c|

% 1) A = LS matrix of linear elasticity
% 2) B x - c>=0   <-->  u n - g <= 0
% where with i in nodes of GammaC
%  c_i= - gap_i
%  B_{i,[2*NE+ n,2*NE+N + n ]}= [ - n_{i,x}, -n_{i,y}]
% here n denotes a node on GammaC

gap=parameters.gap;

normal_node=mesh.normal_node;
normal_node_contact=mesh.normal_node_contact;
normal_edge_contact=mesh.normal_edge_contact;

contC=0;
UN = sparse( length(N_contact), 2 * NE_loc + 2 * N_loc);
ST = sparse( length(E_contact), 2 * NE_loc + 2 * N_loc);
SN = sparse( length(E_contact), 2 * NE_loc + 2 * N_loc);
CN = sparse( length(N_contact), 1);
CT = sparse( length(E_contact), 1);
C=sparse(length(N_contact),1);
contN=0;
for nn=N_contact
    contN=contN+1;
    node=mesh.node(nn,:);
    % if on the boundary and, in particular, on GammaC
        contC=contC+1;
        GammaC(contC)=nn;
        UN(contC,[2*NE_loc + contN, 2*NE_loc + N_loc + contN])=[-normal_node_contact{nn}(1) , -normal_node_contact{nn}(2)] ; 
        
    % C(contC,1)=-gap(node(nn,1), node(nn,2));
        C(contC,1)=-gap_function(node(1), node(2));
end


contCEdge=0;
for ee=1:length(E_contact)
   
        bb=E_contact(ee);
        contCEdge=contCEdge+1;
        n_e_contact=normal_edge_contact{bb};
        phidotn=phi_dot_n(mesh,bb);
        
        % row: multiplier, column: variables
        % n_obstacle (sigma n) <=0
        % n_obstacle= normal of the obstale
        % n = external normal of the mesh
        % phidotn is necessary because has a sign
        % phidotn (n_obstacle_1 C_1 +n_obstacle_2 C_2) <=0
        % -phidotn (n_obstacle_1 C_1 +n_obstacle_2 C_2) >=0
        
        SN(contCEdge, [ee, (ee +NE_loc)] ) = [-phidotn * n_e_contact(1), -phidotn * n_e_contact(2)  ] ;
        CN(contCEdge,1)=0;
        % row: multiplier, column: variables
        % t_obstacle (sigma n) =0 (frictionlesscase)
        % t_obstacle= tangent of the obstale
        % n = external normal of the mesh
        % phidotn is not necessary because we have an equality to zero
        % for the same reason, it does not matter the direction of the
        % tangent
        % phidotn (t_obstacle_1 C_1 +t_obstacle_2 C_2) =0
         %  (n_obstacle_2 C_1 - n_obstacle_1 C_2) =0
        ST(contCEdge, [ee, (ee +NE_loc)] ) = [ n_e_contact(2), -n_e_contact(1)  ] ;        
        CT(contCEdge,1)=0;
        

    
end
% we add the tangent stress constraints into the matrix because they
% areequality constraints. So no need for active set.

STT=ST';

B=[UN;SN];
C=[C;CN;];
BT=B';

end