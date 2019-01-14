function [Patch_Internal_All,Patch_face,Patch_Node]=Arnold_patch(mesh)

L=size(mesh);
L=L(1);

Patch_Boundary_Face=cell(L,1);
Patch_Boundary_Node=cell(L,1);
Patch_Internal_Face=cell(L,1);
Patch_Internal_Node=cell(L,1);


for lev=1:L
N=mesh{lev}.N;
NF=mesh{lev}.NF;
% for each vertex on the mesh of level lev
for n_glob=1:N
    
% compute the number of elements that share that vertex
NT_loc=size(mesh{lev}.N_to_T{n_glob});
NT_loc=NT_loc(2);


% compute the elemFnts that share that vertex
Elems=[];
cont_N=0;
cont_F=0;
cont_all=0;
border_node=[];
border_face=[];
local_face=[];
local_node=[];
% first of all, check if the vertex itself is a vertex-bc
if(mesh{lev}.N_dirichlet(n_glob)==1)
    cont_N=1;
    border_node(cont_N)=n_glob;
end

% then define the elements that share the vertex in Elems
% and find which nodes/vertices are on the boundary of the patch
% all the opposite face/vertices to the node n_glob are boundary
% but also the remaining faces/vertices thare are dirichlet-bc
% Example:
%  4___3
%  |  /|
%  | / |
%  1---|2
% If we have face-dirichlet bc for 12,23,34 and vertex-dirichlet for 1,4.
% If n_glob=1:
% border_face=12,23,34 
% (23,34, because are opposite faces, 12 because is face-dirichlet bc)
% border_node= 1,2,3,4 
%(2,3,4 because are opposite vertices, 1 because is vertec-dirichlet bc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% In 3D a tehtraedron with nodes n1<n2<n3<n4 have the following opposite faces %%%%%%%%%%%%%%%%
%%%% node n1: opposite face n2,n3,n4                                              %%%%%%%%%%%%%%%%
%%%% node n2: opposite face n1,n3,n4                                              %%%%%%%%%%%%%%%%
%%%% node n3: opposite face n1,n2,n4                                              %%%%%%%%%%%%%%%%
%%%% node n4: opposite face n1,n2,n3                                              %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for t_loc=1:NT_loc
Elems(t_loc)=mesh{lev}.N_to_T{n_glob}{t_loc};

T_loc=Elems(t_loc);

% node-dofs of the local t_loc elemFnt
elem=mesh{lev}.elem(T_loc,:);
% face-dofs of the local t_loc elemFnt
elemF=mesh{lev}.elemF(T_loc,:);

F_dirichlet=mesh{lev}.F_dirichlet(elemF);
cont_N=cont_N+3;
cont_F=cont_F+1;

cont_all=cont_all+4;

%%%%%%% qui male
local_face(cont_all-3:cont_all)=[elemF(1);elemF(2);elemF(3);elemF(4)];
local_node(cont_all-3:cont_all)=[elem(1);elem(2);elem(3);elem(4)];
%fixed the vertex n_glob, define the vertices and faces that are on its boundary
        if(elem(1)==n_glob)
            
        border_node(cont_N-2)=elem(2);
        border_node(cont_N-1)=elem(3);
        border_node(cont_N)=elem(4);
        border_face(cont_F)=elemF(1);
        
        otherfaces=[2,3,4];
        for kk=otherfaces
        if(F_dirichlet(kk)==1)
        cont_F=cont_F+1;
        border_face(cont_F)=elemF(kk);    
        end      
        end
        
        elseif(elem(2)==n_glob)

        border_node(cont_N-2)=elem(1);
        border_node(cont_N-1)=elem(3);
        border_node(cont_N)=elem(4);
        border_face(cont_F)=elemF(2);
        
        otherfaces=[1,3,4];
        for kk=otherfaces
        if(F_dirichlet(kk)==1)
        cont_F=cont_F+1;
        border_face(cont_F)=elemF(kk);    
        end        
        end
        
        elseif(elem(3)==n_glob)

        border_node(cont_N-2)=elem(1);
        border_node(cont_N-1)=elem(2);
        border_node(cont_N)=elem(4);
        border_face(cont_F)=elemF(3);
        
        otherfaces=[1,2,4];
        for kk=otherfaces
        if(F_dirichlet(kk)==1)
        cont_F=cont_F+1;
        border_face(cont_F)=elemF(kk);    
        end  
        end
        
        else
            
        border_node(cont_N-2)=elem(1);
        border_node(cont_N-1)=elem(2);
        border_node(cont_N)=elem(3);
        border_face(cont_F)=elemF(4);
        
        otherfaces=[1,2,3];
        for kk=otherfaces
        if(F_dirichlet(kk)==1)
        cont_F=cont_F+1;
        border_face(cont_F)=elemF(kk);    
        end  
        end
        
        end
        
end
% for each level, for each vertex, define the face/vertex border
Patch_Boundary_Face{lev}{n_glob}=unique(border_face);
Patch_Boundary_Node{lev}{n_glob}=unique(border_node);

local_face=unique(local_face);
local_node=unique(local_node);
Patch_face{lev}{n_glob}=local_face;
Patch_Node{lev}{n_glob}=local_node;
% for each level, for each vertex, define the face/vertex internal points
Patch_Internal_Face{lev}{n_glob}=setdiff(local_face,Patch_Boundary_Face{lev}{n_glob});
Patch_Internal_Node{lev}{n_glob}=setdiff(local_node,Patch_Boundary_Node{lev}{n_glob});

vecF=Patch_Internal_Face{lev}{n_glob};
vecN=Patch_Internal_Node{lev}{n_glob};

Patch_Internal_All{lev}{n_glob}=[vecF, ...
                                 vecF + NF,  ...
                                 vecF + 2*NF,  ...
                                 vecN + 3*NF,  ...
                                 vecN + 3*NF + N,  ...
                                 vecN + 3*NF + 2*N];
end





end

end
