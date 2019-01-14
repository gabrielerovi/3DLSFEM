function [A,b,ST,STT,CT,B,BT,C] = ActiveSetBmatrixCvector(mesh,parameters,A,b)
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
L=size(mesh);
L=L(1);
N=mesh{L}.N;
NE=mesh{L}.NE;
N_bc=mesh{L}.N_bc;
E_bc=mesh{L}.E_bc;
node=mesh{L}.node;
edge=mesh{L}.edge;
normal_node=mesh{L}.normal_node;
normal_node_contact=mesh{L}.normal_node_contact;
normal_edge_contact=mesh{L}.normal_edge_contact;

boundary=mesh{L}.boundary;
NB=length(boundary);

[dirichlet,n_and_or_t,bool_bc]= boundary_value_bool(1);

contC=0;
UN = sparse( 0, 2 * NE + 2 * N);
ST = sparse( 0, 2 * NE + 2 * N);
SN = sparse( 0, 2 * NE + 2 * N);
CN = sparse( 0, 1);
CT = sparse( 0, 1);
C=sparse(0,0);
for nn=1:N
    % if on the boundary and, in particular, on GammaC
    if(N_bc(nn)>0 && dirichlet(N_bc(nn),3)==1)
        contC=contC+1;
        GammaC(contC)=nn;
        UN(contC,[2*NE + nn, 2*NE + N + nn])=[-normal_node_contact{nn}(1) , -normal_node_contact{nn}(2)] ;
        
%        C(contC,1)=-gap(node(nn,1), node(nn,2));

        C(contC,1)=-gap_function(node(nn,1), node(nn,2));
    end
    
end

% BT=B';
contCEdge=0;
[dirichlet_E,n_and_or_t,bool_bc]= boundary_value_bool(2);

for bb=1:NE
    % if on the boundary and, in particular, on GammaC
    if(E_bc(bb)>0 && dirichlet_E(E_bc(bb),3)==1)     
        
        contCEdge=contCEdge+1;
        n_e_contact=normal_edge_contact{bb};
        phidotn=phi_dot_n(mesh{L},bb);
        
        % row: multiplier, column: variables
        % n_obstacle (sigma n) <=0
        % n_obstacle= normal of the obstale
        % n = external normal of the mesh
        % phidotn is necessary because has a sign
        % phidotn (n_obstacle_1 C_1 +n_obstacle_2 C_2) <=0
        % -phidotn (n_obstacle_1 C_1 +n_obstacle_2 C_2) >=0
        
        SN(contCEdge, [bb, (bb +NE)] ) = [-phidotn * n_e_contact(1), -phidotn * n_e_contact(2)  ] ;
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
        ST(contCEdge, [bb, (bb +NE)] ) = [ n_e_contact(2), -n_e_contact(1)  ] ;        
        CT(contCEdge,1)=0;
                
    end
    
end


% we add the tangent stress constraints into the matrix because they
% areequality constraints. So no need for active set.

STT=ST';

B=[UN;SN];
C=[C;CN;]
BT=B';

end