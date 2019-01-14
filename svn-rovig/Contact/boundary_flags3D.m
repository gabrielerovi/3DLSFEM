function [N_bc,F_bc,T_bc]=boundary_flags3D(boundary,N,face,F_to_T,face_per_elem)

% boundary is a matrix with all the boundary edges
% N = number of vertices
% edge i
B=length(boundary(:,1));
NF=length(face(:,1));

N_bc=zeros(N,1);
N_cont=zeros(N,1);
F_bc=zeros(NF,1);
T_bc=zeros(B,1);
for bb=1:B
    
    
    % consider the nodes of the edge
    side=sort(boundary(bb,[1:face_per_elem-1]));
    % consider the label of the boundary element
    label=boundary(bb,face_per_elem);
  

    % loop on all the faces to find the one that is on the boundary 
     on_boundary=0;
        for kk=1:NF
        if(face(kk,1:3) == side)
            F_bc( kk )=label;
            T_bc(bb)=F_to_T{kk}{1};
            break;
        end
        end

    % if not, then we impose on both nodes that they are dirichlet-nodal bc
    % otherwise we leave them as zero
     [bool_node,Marker_node,stronger_dirichlet]=is_surface_dirichlet(label,1);
%     [bool,Marker]=is_surface_dirichlet(label,2);
%     [dirichlet,n_and_or_t,bool_bc_node]= boundary_value_bool(1);
    [dirichlet,n_and_or_t,bool_bc]= boundary_value_bool(3);

    for kk=1:face_per_elem-1
        % we have not already found the node or the face is not neumann 
    
    if(N_cont(side(kk))==0 || ~bool_bc(label))

        if(N_cont(side(kk))==1 && stronger_dirichlet==0)
        else
        if(bool_node==1)
          N_bc( (side(kk) ) )= (label); 
          N_cont(side(kk))=1;
        else
          N_bc( (side(kk) ) )= (label);   
        end 
        
        % if we are on a GammaC node, label it as GammaC
        
        if(dirichlet(label,end)==1 &&N_cont(side(kk))==0)
            N_cont(side(kk))=1;
        end
        end
    end
    
  end    
end
 
    T_bc=unique(T_bc);
    
end


