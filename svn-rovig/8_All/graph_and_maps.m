
function [maps,graph]=graph_and_maps(mesh)

[Boundary_Edge,Boundary_Node,EmapLoc2Glob,EmapGlob2Loc,NmapLoc2Glob,NmapGlob2Loc]=local_assembling(mesh);
[Patch_Internal_All,Patch_Edge,Patch_Node,Patch_Boundary_Edge,Patch_Boundary_Node,EmapLoc2Glob,EmapGlob2Loc,NmapLoc2Glob,NmapGlob2Loc]=Arnold_patch(mesh);

L=size(mesh);
L=L(1);
graph=cell(L,1);
for lev=1:L
graph{lev}=graph_neighb(mesh{lev});
end

maps.Boundary_Edge=Boundary_Edge;
maps.Boundary_Node=Boundary_Node;
maps.EmapLoc2Glob=EmapLoc2Glob;
maps.EmapGlob2Loc=EmapGlob2Loc;
maps.NmapLoc2Glob=NmapLoc2Glob;
maps.NmapGlob2Loc=NmapGlob2Loc;
maps.Patch_Edge=Patch_Edge;
maps.Patch_Node=Patch_Node;
maps.Patch_Boundary_Edge=Patch_Boundary_Edge;
maps.Patch_Boundary_Node=Patch_Boundary_Node;
maps.EmapLoc2Glob=EmapLoc2Glob;
maps.EmapGlob2Loc=EmapGlob2Loc;
maps.NmapLoc2Glob=NmapLoc2Glob;
maps.NmapGlob2Loc=NmapGlob2Loc;
maps.Patch_Internal_All=Patch_Internal_All;

end