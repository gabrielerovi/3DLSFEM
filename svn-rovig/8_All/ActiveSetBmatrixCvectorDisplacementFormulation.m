function [B,BT,C] = ActiveSetBmatrixCvectorDisplacementFormulation(mesh,parameters)
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
N_bc=mesh{L}.N_bc;
node=mesh{L}.node;
% consider the normal of the obstacle
normal_node=mesh{L}.normal_node_contact;

[dirichlet,n_and_or_t,bool_bc]= boundary_value_bool(1);

contC=0;
B = sparse( 0, 2 * N);
C=sparse(0,0);
for nn=1:N
    % if on the boundary and, in particular, on GammaC
    if(N_bc(nn)>0 && dirichlet(N_bc(nn),3)==1)
        contC=contC+1;
        GammaC(contC)=nn;
        B(contC,[nn, N + nn])=[-normal_node{nn}(1) , -normal_node{nn}(2)] ;
        C(contC,1)=-gap(node(nn,1), node(nn,2));
        nodeGammaC(nn)=contC;
    else
    nodeGammaC(nn)=0;
    end
    
end

BT=B';






end