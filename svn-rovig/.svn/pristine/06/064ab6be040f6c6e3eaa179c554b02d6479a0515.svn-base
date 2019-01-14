
function [maps,edgegraphinterior,edgegraphgammac]=mapsedge(maps,mesh)

L=length(mesh);
for lev=1:L
    
    edge=mesh{lev}.edge;
    
    for ee=1:length(edge)
        
    % we take all the common faces and the vertex (the last 3 dofs of each
    % non dirichlet node
    maps.Edge_Patch_Internal_All{lev}{ee}=sort(intersect(maps.Patch_Internal_All{lev}{edge(ee,1)},maps.Patch_Internal_All{lev}{edge(ee,2)}));
    for jj=1:2
    if(mesh{lev}.N_dirichlet(edge(ee,jj))==0)
    maps.Edge_Patch_Internal_All{lev}{ee}=unique([maps.Edge_Patch_Internal_All{lev}{ee},maps.Patch_Internal_All{lev}{edge(ee,jj)}(end-2:end)]);
    end
    end

    end
    
    [edgegraphinterior{lev},edgegraphgammac{lev}]=graph_neighb_interior_and_gammac(mesh{lev});
end
end