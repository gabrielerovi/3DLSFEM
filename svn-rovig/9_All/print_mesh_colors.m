function print_mesh_colors(mesh,nodes_colors)
L=length(mesh);
for lev=1:L
figure
N=mesh{lev}.N;

radius(1)=(max(mesh{lev}.node(:,1))-min(mesh{lev}.node(:,1)))/(N);
radius(2)=(max(mesh{lev}.node(:,2))-min(mesh{lev}.node(:,2)))/(N);
radius=min(radius);
number_of_color=max(nodes_colors{lev});
for ii=1:mesh{lev}.N
    
center = [mesh{lev}.node(ii,1),mesh{lev}.node(ii,2)];
if(nodes_colors{lev}(ii)==1)
    local_color=[1,0,0];
elseif(nodes_colors{lev}(ii)==2)
    local_color=[0,1,0];
elseif(nodes_colors{lev}(ii)==3)
    local_color=[0,0,1];
elseif(nodes_colors{lev}(ii)==4)
    local_color=[0.5,0.2,0];
elseif(nodes_colors{lev}(ii)==5)
    local_color=[0,0.5,0.2];
elseif(nodes_colors{lev}(ii)==6)
    local_color=[0.2,0,0.5];
elseif(nodes_colors{lev}(ii)==7)
    local_color=[0.3,0.3,0];    
elseif(nodes_colors{lev}(ii)==7)
    local_color=[0,0.3,0.3];  
elseif(nodes_colors{lev}(ii)==7)
    local_color=[0.3,0.,0.3];  
else    
    local_color=[value,value,value];
end
viscircles(center,radius,'Color',local_color);
end

for ii=1:mesh{lev}.NE
side=mesh{lev}.edge(ii,:);
node1=mesh{lev}.node(side(1),:);
node2=mesh{lev}.node(side(2),:);


coeff1=1;%NDtoRT{L}(ii,mesh{L}.edge(ii,1));
coeff2=1;%NDtoRT{L}(ii,mesh{L}.edge(ii,2));

x = 0.5 * (node1(1)+node2(1));
y = 0.5 * (node1(2)+node2(2));
txt.Color='b';  

hold on
arrow = 0.8*(node2-node1);                         % Difference

quiver(node1(1),node1(2),arrow(1),arrow(2),0)


end
end


end