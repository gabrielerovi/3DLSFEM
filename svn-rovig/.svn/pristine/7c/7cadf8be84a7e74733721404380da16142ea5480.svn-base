
function [maps,graph]=graph_and_maps(mesh)

[Patch_Internal_All,Patch_Face,Patch_Node]=Arnold_patch(mesh);

L=size(mesh);
L=L(1);
graph=cell(L,1);
for lev=1:L
graph{lev}=graph_neighb(mesh{lev});
% graph{lev}=graph_neighb_interior_and_gammac(mesh{lev});
end

maps.Patch_Face=Patch_Face;
maps.Patch_Node=Patch_Node;
maps.Patch_Internal_All=Patch_Internal_All;

end