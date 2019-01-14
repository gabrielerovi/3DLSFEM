function [nodes_colors,colors_cardinality]=UTOPIAcoloring(mesh,maps,coloring_sequential)


L=length(mesh);
for lev=1:L
N=mesh{lev}.N;
nodes_colors{lev}=zeros(N,1);
colors=[];
num_colors=0;
colors_cardinality=cell(1,1);
colors_cardinality{1}=zeros(1,1);

nodes_colors{lev}(1)=1;
colors_cardinality{lev}(1)=1; 
if(coloring_sequential)
for nn=2:N
    
    
    %%%% patch_colors contains all the colors not already used in the patch
    %%%% setdiff(1:num_colors,patch_colors) contains the colors that can be
    %%%% used, among which we take the one with minimum colors_cardinality
    %%%% if patch_colors=[], then all the colors have been used and we must
    %%%% add a new colorc
    patch_colors=setdiff(1:num_colors,nodes_colors{lev}(maps.Patch_Node{lev}{nn}));
    % if there is no color left, then add an extra color
    if( isempty(patch_colors))
    num_colors=num_colors+1;
    colors_cardinality{lev}(num_colors)=1;
    nodes_colors{lev}(nn)=num_colors;  
    % otherwise consider the color with minimum cardinality
    else
    [min_val,min_pos]=min(colors_cardinality{lev});
    nodes_colors{lev}(nn)=min_pos;   
    colors_cardinality{lev}(min_pos)=colors_cardinality{lev}(min_pos)+1;
    end
    
end
else
cont=1;
num_colors=1;
already_colored=zeros(N,1);
already_colored(1)=1;
patch=maps.Patch_Node{lev}{1};
tmp=setdiff(patch,1);
neighb_nn=tmp(1);
nn=neighb_nn(1);
while(cont<N) 
    patch_colors=setdiff(1:num_colors,nodes_colors{lev}(maps.Patch_Node{lev}{nn}));
    if( isempty(patch_colors))
    num_colors=num_colors+1;
    colors_cardinality{lev}(num_colors)=1;
    nodes_colors{lev}(nn)=num_colors;  
    % otherwise consider the color with minimum cardinality
    else
    colors_cardinality_loc=colors_cardinality{lev}(patch_colors);

    [min_val,min_pos]=min(colors_cardinality_loc);
    nodes_colors{lev}(nn)=patch_colors(min_pos);   
    colors_cardinality{lev}(min_pos)=colors_cardinality{lev}(min_pos)+1;
    end    
    already_colored(nn)=1;
    patch=maps.Patch_Node{lev}{nn};
    neighb_nn=patch(find(already_colored(patch)==0));
    if(isempty(neighb_nn))
        nn=min(find(already_colored==0));
    else
        nn=neighb_nn(1);
    end
    cont=cont+1;
    
end   
end


end
end
