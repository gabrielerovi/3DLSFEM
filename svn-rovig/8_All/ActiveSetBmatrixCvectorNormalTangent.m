function [A,b,ST,STT,CT,B,BT,C] = ActiveSetBmatrixCvectorNormalTangent(mesh,parameters,R,RT,A,b)
% we solve for min H, with H=0.5 x' A x - x' f - lambda (B x -c)
% structure of the problem
% |A -B'| |x     |= |f|
% |B  0 | |lambda|  |c|

% 1) A = LS matrix of LS linear elasticity
% 2) B x - c>=0   <-->  u n - g <= 0 AND n_obst' sigma n <=0
% here we consider on gammaC dofs that have been transformed into normal
% and tangent components. Therefor we just have u_n<=g, s_n <=0. In the
% matrix B we will have only zeros and ones.

gap=parameters.gap;
L=size(mesh);
L=L(1);
N=mesh{L}.N;
NE=mesh{L}.NE;
node=mesh{L}.node;

contC=0;
UN = sparse( 0, 2 * NE + 2 * N);
ST = sparse( 0, 2 * NE + 2 * N);
SN = sparse( 0, 2 * NE + 2 * N);
CN = sparse( 0, 1);
CT = sparse( 0, 1);
C=sparse(0,0);

% the first component is normal, the second one is tangent
for nn=mesh{L}.N_contact
    % if on the boundary and, in particular, on GammaC
        contC=contC+1;
        GammaC(contC)=nn;
        UN(contC,2*NE + nn)=-1 ;
        
%        C(contC,1)=-gap(node(nn,1), node(nn,2));

        C(contC,1)=-gap_function(node(nn,1), node(nn,2));
end

contCEdge=0;
for bb=mesh{L}.E_contact
    % if on the boundary and, in particular, on GammaC        
        contCEdge=contCEdge+1;
        % row: multiplier, column: variables
        % n_obstacle (sigma n) <=0
        % n_obstacle= normal of the obstale
        % n = external normal of the mesh
        % phidotn is necessary because has a sign
        % phidotn (n_obstacle_1 C_1 +n_obstacle_2 C_2) <=0
        % -phidotn (n_obstacle_1 C_1 +n_obstacle_2 C_2) >=0
        phidotn=phi_dot_n(mesh{L},bb);
        SN(contCEdge, bb  ) = -1 * phidotn;
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
        ST(contCEdge, bb +NE ) =1 ;        
        CT(contCEdge,1)=0;
    
end

% we add the tangent stress constraints into the matrix because they
% areequality constraints. So no need for active set.

STT=ST';

B=[UN;SN];
C=[C;CN;]
BT=B';

A=R*A*RT;
b=R*b;


end