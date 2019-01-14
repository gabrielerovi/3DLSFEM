function [N_bc,E_bc,F_bc,T_bc]=boundary_flags3D(boundary, N, edge, face,F_to_T)

% boundary is a matrix with all the boundary faces
% N = number of vertices
% edge i
B=length(boundary(:,1));
NE=length(edge(:,1));
NF=length(face(:,1));

N_bc=zeros(N,1);
E_bc=zeros(NE,1);
F_bc=zeros(NF,1);
T_bc=zeros(NF,1);
vertex_per_face=3;
edge_per_face=3;
for bb=1:B
    
    
    % consider the nodes of the faces
    vertex=sort(boundary(bb,[1:vertex_per_face]));
    edges(1,[1,2])=[vertex(1), vertex(2)];
    edges(2,[1,2])=[vertex(2), vertex(3)];
    edges(3,[1,2])=[vertex(1), vertex(3)];    
    faces=[vertex(1), vertex(2),vertex(3)];
    
    % consider the label of the boundary element
    label=boundary(bb,vertex_per_face+1);
  
    % loop on the nodes of the boundary element and impose label
    for kk=1:vertex_per_face
        if( N_bc(  (vertex(kk) ) )==0)
            N_bc( (vertex(kk) ) )=label;
        end
    end
    % loop on all the edges to find the one that is on the boundary 
    for ee=1:edge_per_face
        for kk=1:NE
        if(edge(kk,[1,2]) == edges(ee,:))
            E_bc( kk )=label;
            break;
        end
        end
    end
    % loop on all the edges to find the one that is on the boundary 
        for kk=1:NF
        if(face(kk,[1,2,3]) == faces)
            F_bc( kk )=label;
            T_bc(bb)=F_to_T{kk}{1};
            break;
        end
        end
end
 
    
end


