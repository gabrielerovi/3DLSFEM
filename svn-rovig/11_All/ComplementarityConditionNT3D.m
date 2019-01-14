function [ComplementarityMatrix,Complementarityb] =ComplementarityConditionNT3D(L,mesh,parameters)
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
NF=mesh{L}.NF;
N=mesh{L}.N;
dim=parameters.dim;
ComplementarityMatrix=sparse( (N+NF) * dim , (N+NF)*dim );
Complementarityb=sparse(  (N+NF) * dim ,1);
for ff=mesh{L}.F_contact
                
            % find the node dofs of the face
            face=mesh{L}.face(ff,:);
            % find the coordinates of the node
            nodes_coord_of_face=mesh{L}.node(face,:);
            Area=AreaTriangle(nodes_coord_of_face); 
            

            temporary= parameters.C_contact * 1.0/3 * Area;
            
            for kk=1:3
            g(kk) = gap_function(nodes_coord_of_face(kk,1),nodes_coord_of_face(kk,2),nodes_coord_of_face(kk,3)) ; 
            end
            
            Complementarityb(ff) = temporary * sum(g);

            for ii=1:3
                ComplementarityMatrix(ff,3*NF+face(ii))=temporary;
                ComplementarityMatrix(3*NF+face(ii),ff)=temporary;
            end
            
end
    
end








