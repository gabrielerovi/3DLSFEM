function [Y,W]=dirichlet_bc_vector(X,type_of_dof,mesh)

% X contains the labels

N=length(X);
Y=zeros(N,2);
W=[];
[dirichlet,n_and_or_t,bool]= boundary_value_bool(type_of_dof);
% dirichlet_node=[0.05 0.0 0.0 0.0 0 0 ];
% dirichlet_edge=[0. 0 0 0 0 0];
cont=0;
% dof=node
if(type_of_dof==1)    
    for ii=1:N
        [bool_loc,Marker_loc]=is_surface_dirichlet(X(ii),type_of_dof);
        if(bool_loc>0)
        cont=cont+1;
        Y(ii,:) = dirichlet(X(ii),[1,2]);
        W(cont,:) = dirichlet(X(ii),[1,2])
        end
    end
else if(type_of_dof==2)
        kk(1)=3;kk(2)=1;kk(3)=2;
    for ii=1:N
        [bool_loc,Marker_loc]=is_surface_dirichlet(X(ii),type_of_dof);
        if(bool_loc>0)
        % and edge on the boundary in 2D belongs to only 1 element
        E_to_T=cell2mat(mesh.E_to_T{ii});
        elemE=mesh.elemE(E_to_T,:);
        nodeE=mesh.node(mesh.elem(E_to_T,:),:);  
        for mm=1:3
            if(elemE(mm)==ii)
                cc=setdiff(kk,kk(mm));
                break;
            end
        end
        midpoint=sum(nodeE(cc,:))/length(cc);
        opposite_to_midpoint=midpoint-nodeE(kk(mm),:);
        side=mesh.node(mesh.edge(ii,:),:);
        RT_normal(1)=side(2,2)-side(1,2);
        RT_normal(2)=-side(2,1)+side(1,1);        
        if(opposite_to_midpoint*RT_normal'>0)
            sign=1;
        else
            sign=-1;
        end
        coeff=sign*norm(side(1,:)-side(2,:));
        cont=cont+1;
        Y(ii,:) = coeff*dirichlet(X(ii),[1,2])
        W(cont,:) = coeff*dirichlet(X(ii),[1,2])
        end
    end    

end

end