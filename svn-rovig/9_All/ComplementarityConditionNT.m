function [ComplementarityMatrix,Complementarityb] =ComplementarityConditionNT(L,mesh,parameters)
%%%% We want to discretize < sigma_n, u_n -g  >_GammaC
% u_n is discretized linearly, so on an edge ee of point n1 and n2
% u_n=t u_n1 + (1-t) u_n2
% sigma_n is constant on each edge. It is derived from the following
% transformation: [sigma_n, sigma_t] = O [f_x,f_y], where O is the matrix
% that transform from x-y to n-t coordinates.
% Additionally we have to convert:
%[f_x,f_y]=sigma n =[phi S1, phi S2 ]' [n1, n2]=(phi n) [S1,S2 ]'
% Therefore: [s_n,s_t]= O (phi n) [S1,S2 ]' = H [S1, S2]
% Since all dofs on GammaC are transformed with the HouseHolder
% transformation H, then our dofs on GammaC are only [s_n,s_t].
% Also the firs component is the normal one and it is constant on the whole
% edge. Therefore the contribution from s_n on an edge ee is simply 1.

% int_E simga_n (u_n-g)dE=int_E phi_sn sigma_n ( u_n1 t + (1-t) u_n2 - g1 t- (1-t) g2)
% where phi_sn=1, so:
% V1 = sigma_n u_n1  int_E    t  dE = sigma_n u_n1 0.5 L
% V2 = sigma_n u_n2  int_E (1-t) dE = sigma_n u_n2 0.5 L
% F1 = sigma_n g1    int_E g1 t + g2 (1-t) dE = sigma_n g1 0.5 L (g1+g2)
NE=mesh{L}.NE;
N=mesh{L}.N;

ComplementarityMatrix=sparse(NE*2+N*2,NE*2+N*2);
Complementarityb=sparse(NE*2+N*2,1);
for ee=mesh{L}.E_contact
            % find the node dofs of the edge
            edge=mesh{L}.edge(ee,:);
            % find the coordinates of the node
            point=mesh{L}.node(edge,:);
            
            
            % the normals we are going to use are the ones of the obstacle
%             normal_edge=mesh{L}.normal_edge_contact{ee};
%             normal_node(1,:)=mesh{L}.normal_node_contact{edge(1)};
%             normal_node(2,:)=mesh{L}.normal_node_contact{edge(2)};
%                         
            halfLside= parameters.C_contact * 0.5 * norm( point(1,1:2)-point(2,1:2));
            
            g1 = gap_function(point(1,1),point(1,2)) ; 
            g2 = gap_function(point(2,1),point(2,2)) ;
            
            Complementarityb(ee) = halfLside * ( g1+g2);

            
            for ii=1:2
                ComplementarityMatrix(ee,2*NE+edge(ii))=halfLside;
                ComplementarityMatrix(2*NE+edge(ii),ee)=halfLside;
%             A{L,1,3}(ee,edge(ii))=A{L,1,3}(ee,edge(ii))+halfLside;
%             A{L,3,1}(edge(ii),ee)=A{L,3,1}(edge(ii),ee)+halfLside;
            end
end
    
end








