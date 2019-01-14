function [Patch_Internal_All,Patch_Edge,Patch_Node,Patch_Boundary_Edge,Patch_Boundary_Node,EmapLoc2Glob,EmapGlob2Loc,NmapLoc2Glob,NmapGlob2Loc]=Arnold_patch(mesh)

L=size(mesh);
L=L(1);

Patch_Boundary_Edge=cell(L,1);
Patch_Boundary_Node=cell(L,1);
Patch_Internal_Edge=cell(L,1);
Patch_Internal_Node=cell(L,1);
EmapLoc2Glob=cell(L,1);
EmapGlob2Loc=cell(L,1);
NmapLoc2Glob=cell(L,1);
NmapGlob2Loc=cell(L,1);

for lev=1:L
N=mesh{lev}.N;
NE=mesh{lev}.NE;
% for each vertex on the mesh of level lev
for n_glob=1:N
    
% compute the number of elements that share that vertex
NT_loc=size(mesh{lev}.N_to_T{n_glob});
NT_loc=NT_loc(2);


% compute the elements that share that vertex
Elems=[];
cont_N=0;
cont_E=0;
cont_all=0;
border_node=[];
border_edge=[];
local_edge=[];
local_node=[];
% first of all, check if the vertex itself is a vertex-bc
if(mesh{lev}.N_dirichlet(n_glob)==1)
    cont_N=1;
    border_node(cont_N)=n_glob;
end

% then define the elements that share the vertex in Elems
% and find which nodes/vertices are on the boundary of the patch
% all the opposite edge/vertices to the node n_glob are boundary
% but also the remaining edges/vertices thare are dirichlet-bc
% Example:
%  4___3
%  |  /|
%  | / |
%  1---|2
% If we have edge-dirichlet bc for 12,23,34 and vertex-dirichlet for 1,4.
% If n_glob=1:
% border_edge=12,23,34 
% (23,34, because are opposite edges, 12 because is edge-dirichlet bc)
% border_node= 1,2,3,4 
%(2,3,4 because are opposite vertices, 1 because is vertec-dirichlet bc)
for t_loc=1:NT_loc
Elems(t_loc)=mesh{lev}.N_to_T{n_glob}{t_loc};

T_loc=Elems(t_loc);

% node-dofs of the local t_loc element
elem=mesh{lev}.elem(T_loc,:);
% edge-dofs of the local t_loc element
elemE=mesh{lev}.elemE(T_loc,:);

E_dirichlet=mesh{lev}.E_dirichlet(elemE);
cont_N=cont_N+2;
cont_E=cont_E+1;

cont_all=cont_all+3;
local_edge(cont_all-2:cont_all)=[elemE(1);elemE(2);elemE(3)];
local_node(cont_all-2:cont_all)=[elem(1);elem(2);elem(3)];
%fixed the vertex n_glob, define the vertices and edges that are on its boundary
        if(elem(1)==n_glob)
        border_node(cont_N-1)=elem(2);
        border_node(cont_N)=elem(3);
        border_edge(cont_E)=elemE(2);
        if(E_dirichlet(1)==1)
        cont_E=cont_E+1;
        border_edge(cont_E)=elemE(1);    
        end
        if(E_dirichlet(3)>0)
        cont_E=cont_E+1;
        border_edge(cont_E)=elemE(3); 
        end
        elseif(elem(2)==n_glob)
        border_node(cont_N-1)=elem(1);
        border_node(cont_N)=elem(3);
        border_edge(cont_E)=elemE(3);   
        if(E_dirichlet(1)==1)
        cont_E=cont_E+1;
        border_edge(cont_E)=elemE(1);    
        end
        if(E_dirichlet(2)==1)
        cont_E=cont_E+1;
        border_edge(cont_E)=elemE(2); 
        end        
        else
        border_node(cont_N-1)=elem(1);
        border_node(cont_N)=elem(2);
        border_edge(cont_E)=elemE(1);  
        if(E_dirichlet(2)==1)
        cont_E=cont_E+1;
        border_edge(cont_E)=elemE(2);    
        end
        if(E_dirichlet(3)==1)
        cont_E=cont_E+1;
        border_edge(cont_E)=elemE(3); 
        end        
        end
        
end
% for each level, for each vertex, define the edge/vertex border
Patch_Boundary_Edge{lev}{n_glob}=unique(border_edge);
Patch_Boundary_Node{lev}{n_glob}=unique(border_node);

local_edge=unique(local_edge);
local_node=unique(local_node);
Patch_Edge{lev}{n_glob}=local_edge;
Patch_Node{lev}{n_glob}=local_node;
% for each level, for each vertex, define the edge/vertex internal points
Patch_Internal_Edge{lev}{n_glob}=setdiff(local_edge,Patch_Boundary_Edge{lev}{n_glob});
Patch_Internal_Node{lev}{n_glob}=setdiff(local_node,Patch_Boundary_Node{lev}{n_glob});

vecE=Patch_Internal_Edge{lev}{n_glob};
vecN=Patch_Internal_Node{lev}{n_glob};

Patch_Internal_All{lev}{n_glob}=[vecE, vecE+NE, vecN+2*NE, vecN+2*NE + N];

% for each level, for each vertex, define the map:
% EmapGlob2Loc: global edge/vertex -> global edge/vertex
% EmapLoc2Glob: local edge/vertex -> global edge/vertex
E_map=unique(mesh{lev}.elemE(Elems,:));
N_map=unique(mesh{lev}.elem(Elems,:));
EmapGlob2Loc{lev}{n_glob} = containers.Map(E_map,1:length(E_map));
EmapLoc2Glob{lev}{n_glob} = containers.Map(1:length(E_map),E_map);
NmapGlob2Loc{lev}{n_glob} = containers.Map(N_map,1:length(N_map));
NmapLoc2Glob{lev}{n_glob} = containers.Map(1:length(N_map),N_map);


end





end

end
