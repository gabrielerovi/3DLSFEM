function [Boundary_Edge,Boundary_Node,EmapLoc2Glob,EmapGlob2Loc,NmapLoc2Glob,NmapGlob2Loc]=local_assembling(mesh)

L=size(mesh);
L=L(1);

Boundary_Edge=cell(L,1);
Boundary_Node=cell(L,1);
EmapLoc2Glob=cell(L,1);
EmapGlob2Loc=cell(L,1);
NmapLoc2Glob=cell(L,1);
NmapGlob2Loc=cell(L,1);

for lev=1:L
N=mesh{lev}.N;

% for each vertex on the mesh of level lev
for n_glob=1:N
    
% compute the number of elements that share that vertex
NT_loc=size(mesh{lev}.N_to_T{n_glob});
NT_loc=NT_loc(2);


% compute the elements that share that vertex
Elems=[];
cont_N=0;
cont_E=0;
border_node=[];
border_edge=[];
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

elem=mesh{lev}.elem(T_loc,:);
elemE=mesh{lev}.elemE(T_loc,:);
E_dirichlet=mesh{lev}.E_dirichlet(elemE);
cont_N=cont_N+2;
cont_E=cont_E+1;
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
Boundary_Edge{lev}{n_glob}=unique(border_edge);
Boundary_Node{lev}{n_glob}=unique(border_node);

% for each level, for each vertex, define the map:
% EmapGlob2Loc: global edge/vertex -> global edge/vertex
% EmapLoc2Glob: local edge/vertex -> global edge/vertex
E_map=unique(mesh{lev}.elemE(Elems,:));
N_map=unique(mesh{lev}.elem(Elems,:));
EmapGlob2Loc{lev}{n_glob} = containers.Map(E_map,1:length(E_map));
EmapLoc2Glob{lev}{n_glob} = containers.Map(1:length(E_map),E_map);
NmapGlob2Loc{lev}{n_glob} = containers.Map(N_map,1:length(N_map));
NmapLoc2Glob{lev}{n_glob} = containers.Map(1:length(N_map),N_map);


% compute the local system
% if the vertex corresponding to n_glob is internal: 
% # local vertex= 1 + NT_loc,  # local edge = 2 NT_loc
% if the vertex corresponding to n_glob is on the boundary: 
% # local vertex= 2 + NT_loc,  # local edge = 2 NT_loc+1

if(mesh{lev}.N_bc(n_glob)>0)
    NE_loc= 2 * NT_loc + 1;
    N_loc = NT_loc + 2;
else
    NE_loc= 2 * NT_loc;
    N_loc = NT_loc + 1;
end



end

end

end
