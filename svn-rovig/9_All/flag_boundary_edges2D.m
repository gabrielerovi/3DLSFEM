function [flag_for_bc_edge,flag_for_each_edge] =flag_boundary_edges2D  (boundary,edge)

NE=length(edge(:,1));
BE=length(boundary(:,1));

flag_for_bc_edge=zeros(BE,1);
flag_for_each_edge=zeros(NE,1);
cont=0;

for ee=1:NE
    
    for bb=1:BE
        
        if(sort(boundary(bb,[1,2])) == edge (ee,:))
            cont=cont+1;
            flag_for_bc_edge(cont)=boundary(bb,3);
            flag_for_each_edge(ee)=boundary(bb,3);
        end
    end
    
end

end