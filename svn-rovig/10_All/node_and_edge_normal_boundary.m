function [mesh]=node_and_edge_normal_boundary(mesh,parameters)
L=length(mesh);
kk=[3;1;2];
for lev=1:L
grid=mesh{lev};
node_per_elem=grid.node_per_elem;
edge_per_elem=grid.edge_per_elem;
node=grid.node;
edge=grid.edge;
NT=grid.NT;
N=grid.N;
NE=grid.NE;
E_bc=grid.E_bc;

[dirichlet_E,n_and_or_t_E,bool_bc_E]= boundary_value_bool(2);


for nn=1:N
%     normal_node{lev,nn}=[0;0];
%     normal_node_patch{lev,nn}=[];
%     normal_node_contact{lev,nn}=[0;0];
end
for ee=1:NE
    normal_edge{lev,ee}=[];
%     normal_edge_contact{lev,ee}=[0;0];
end

for tt=1:NT
elemE=grid.elemE(tt,:);
elem=grid.elem(tt,:);


for ee=1:edge_per_elem
    eetot=elemE(ee);
    if(E_bc(eetot)>0)
        vertices=edge(eetot,:);
        opposite_vertex=setdiff(elem,edge(eetot,:));
        side=node(vertices,:);
        normal_edge{lev,eetot}=[ side(2,2)-side(1,2);
                         -side(2,1)+side(1,1)];
        normal_edge{lev,eetot}=normal_edge{lev,eetot}/norm(normal_edge{lev,eetot});
        
        tmp=side(1,[1,2])-node(opposite_vertex,[1,2]);
        if(tmp*normal_edge{lev,eetot}>0)
        sign=1;
        else
        sign=-1;
        end
        normal_edge{lev,eetot}=sign*normal_edge{lev,eetot};
%         normal_node{lev,vertices(1)}=normal_node{lev,vertices(1)}+normal_edge{lev,eetot};
%         normal_node{lev,vertices(2)}=normal_node{lev,vertices(2)}+normal_edge{lev,eetot};
%         normal_node_patch{lev,vertices(1)}=[normal_node_patch{lev,vertices(1)},normal_edge{lev,eetot};]    ;
%         normal_node_patch{lev,vertices(2)}=[normal_node_patch{lev,vertices(2)},normal_edge{lev,eetot};]    ;
        
        % if we are on GammaC, then we sum up to the previous normal, the
        % new one
%         if(dirichlet_E(E_bc(eetot),end)==1 && parameters.body_normal_bool==true)
%         normal_node_contact{lev,vertices(1)}=normal_node_contact{lev,vertices(1)}+normal_edge{lev,eetot};
%         normal_node_contact{lev,vertices(2)}=normal_node_contact{lev,vertices(2)}+normal_edge{lev,eetot};
%         normal_edge_contact{lev,eetot}= normal_edge{lev,eetot};
%         elseif(dirichlet_E(E_bc(eetot),end)==1 && parameters.body_normal_bool==false)
%             
%         normal_node_contact{lev,vertices(1)}=[ parameters.body_normal_x( side(1,1), side(1,2)); parameters.body_normal_y( side(1,1), side(1,2))] ;
%         normal_node_contact{lev,vertices(2)}=[ parameters.body_normal_x( side(2,1), side(2,2)); parameters.body_normal_y( side(2,1), side(2,2))] ;
%         normal_node_contact{lev,vertices(1)}=normal_node_contact{lev,vertices(1)}/norm(normal_node_contact{lev,vertices(1)});
%         normal_node_contact{lev,vertices(2)}=normal_node_contact{lev,vertices(2)}/norm(normal_node_contact{lev,vertices(2)});
%         midpointx=0.5* (side(1,1)+side(2,1) );
%         midpointy=0.5* (side(1,2)+side(2,2) );
%         normal_edge_contact{lev,eetot}=[ parameters.body_normal_x( midpointx, midpointy); parameters.body_normal_y( midpointx, midpointy)] ;
%         normal_edge_contact{lev,eetot}=normal_edge_contact{lev,eetot}/norm(normal_edge_contact{lev,eetot});
%         end
    end
end
end

% for nn=1:N
%    if(isempty(normal_node{lev,nn})==false)
%        normal_node{lev,nn}=normal_node{lev,nn}/norm(normal_node{lev,nn});
%        normal_node_contact{lev,nn}=normal_node_contact{lev,nn}/norm(normal_node_contact{lev,nn});
%    end
% end
end


for lev=1:L
    for nn=1:mesh{lev}.N
%     mesh{lev}.normal_node{nn}=normal_node{lev,nn};
%     mesh{lev}.normal_node_contact{nn}=normal_node_contact{lev,nn};
%     mesh{lev}.normal_node_patch{nn}=normal_node_patch{lev,nn};
    end
    for ee=1:mesh{lev}.NE
    mesh{lev}.normal_edge{ee}=normal_edge{lev,ee};
%      mesh{lev}.normal_edge_contact{ee}=normal_edge_contact{lev,ee};
    end
end
end