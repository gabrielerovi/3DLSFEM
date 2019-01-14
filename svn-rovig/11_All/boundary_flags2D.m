function [N_bc,E_bc,T_bc]=boundary_flags2D(boundary, N, edge,E_to_T)

% boundary is a matrix with all the boundary edges
% N = number of vertices
% edge i
B=length(boundary(:,1));
NE=length(edge(:,1));

N_bc=zeros(N,1);
N_cont=zeros(N,1);
E_bc=zeros(NE,1);
T_bc=zeros(B,1);
for bb=1:B
    
    
    % consider the nodes of the edge
    side=sort(boundary(bb,[1:2]));
    % consider the label of the boundary element
    label=boundary(bb,2+1);
  

    % loop on all the edges to find the one that is on the boundary 
     on_boundary=0;
        for kk=1:NE
        if(edge(kk,[1,2]) == side)
            E_bc( kk )=label;
            T_bc(bb)=E_to_T{kk}{1};
            break;
        end
        end

    % if not, then we impose on both nodes that they are dirichlet-nodal bc
    % otherwise we leave them as zero
    [bool_node,Marker_node]=is_surface_dirichlet(label,1);
    [bool,Marker]=is_surface_dirichlet(label,2);
    [dirichlet,n_and_or_t,bool_bc_node]= boundary_value_bool(1);
    [dirichlet,n_and_or_t,bool_bc]= boundary_value_bool(2);

    for kk=1:2
    if(N_cont(side(kk))==0)
        if(bool_node==1)
          N_bc( (side(kk) ) )= (label); 
          N_cont(side(kk))=1;
        else
          N_bc( (side(kk) ) )= (label);   
        end 
        
        % if we are on a GammaC node, label it as GammaC
        
        if(dirichlet(label,3)==1)
            N_cont(side(kk))=1;
        end
    end
    end
    
 

    
    
    
    
%      if(bool1==1 && bool==0)
%         for kk=1:2
%         N_bc( (side(kk) ) )= label; 
%         end
%       % check if the edge is a dirichlet RT-bc
%       N_cont(side)=N_cont(side)+[1;1];
%      end
    
     
     
%     for kk=1:2
%     % this is the second time we find the node on the boundary
%     % and this time it is a dirichlet boundary
%     % so we overwrite the old value
%     if(N_cont(side(kk))==1)
%          N_bc( (side(kk) ) )= label;         
%     elseif(bool==0)
%          N_bc( (side(kk) ) )= label; 
%     end    
%     end
    

    % loop on the nodes of the boundary element and impose label
    % we impose the biggest label. This means that Dirichlet BC label have
    % to be greater than Neumann BC label 
%     for kk=1:2
%             N_bc( (side(kk) ) )=max(N_bc( (side(kk) ) ),label);
%     end
    
    
end
 
    T_bc=unique(T_bc);
    
end


