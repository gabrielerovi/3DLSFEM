
function [maps,graph]=graph_and_maps(mesh)

[Patch_Internal_All,Patch_Edge,Patch_Node]=Arnold_patch(mesh);

L=size(mesh);
L=L(1);
graph=cell(L,1);
for lev=1:L
graph{lev}=graph_neighb(mesh{lev});
end

maps.Patch_Edge=Patch_Edge;
maps.Patch_Node=Patch_Node;
maps.Patch_Internal_All=Patch_Internal_All;

end