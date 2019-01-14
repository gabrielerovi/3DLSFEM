function coeff = RT_dirichlet_coeff(Edge, mesh)

% X contains the labels
        % Edge contains the dof-edge (usefull after the renumbering used to
        % put the bc above in the system)
        type_of_dof=2;
        kk(1)=3;kk(2)=1;kk(3)=2;
        
%         [bool_loc,Marker_loc]=is_surface_dirichlet(X,type_of_dof);
%         if(bool_loc>0)
        % and edge on the boundary in 2D belongs to only 1 element
        %consider the element relative to this edge
        E_to_T=cell2mat(mesh.E_to_T{Edge});
        % consider the edges of the element
        elemE=mesh.elemE(E_to_T,:);
        % consider the nodes of the element
        nodeE=mesh.node(mesh.elem(E_to_T,:),[1,2]);  
        for mm=1:3
            if(elemE(mm)==Edge)
                % cc are the other two edges different from the actual Edge
                cc=setdiff(kk,kk(mm));
                break;
            end
        end
        % midpoint contains the midpoints
        midpoint=sum(nodeE(cc,:))/length(cc);
        % consider the vector from its opposite vertex to itself  
        opposite_to_midpoint=midpoint-nodeE(kk(mm),:);
        % consider the vertices of the actual Edge
        side=mesh.node(mesh.edge(Edge,:),:);
        % define the RT normal
        RT_normal(1)=side(2,2)-side(1,2);
        RT_normal(2)=-side(2,1)+side(1,1);  
        % define the sign
        if(opposite_to_midpoint*RT_normal'>0)
            sign=1;
        else
            sign=-1;
        end
        % define the coefficient
        coeff=sign/norm(side(1,:)-side(2,:));
   
end
      