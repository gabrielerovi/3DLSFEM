function print_mesh_L(mesh,uwant2contour,contour_dof)
L=length(mesh);
for lev=1:L
figure
for ii=1:mesh{lev}.N
    
x = mesh{lev}.node(ii,1);
y = mesh{lev}.node(ii,2);
txt = text(x,y,num2str(ii));
txt.Color='r';  
end

for ii=1:mesh{lev}.NE
side=mesh{lev}.edge(ii,:);
node1=mesh{lev}.node(side(1),:);
node2=mesh{lev}.node(side(2),:);


coeff1=1;%NDtoRT{L}(ii,mesh{L}.edge(ii,1));
coeff2=1;%NDtoRT{L}(ii,mesh{L}.edge(ii,2));

x = 0.5 * (node1(1)+node2(1));
y = 0.5 * (node1(2)+node2(2));
txt = text(x,y,num2str(ii)) ;
txt.Color='b';  

hold on
arrow = 0.8*(node2-node1);                         % Difference

quiver(node1(1),node1(2),arrow(1),arrow(2),0)


end

if(uwant2contour{lev})
cdN=contour_dof{lev}(find(contour_dof{lev}>2*mesh{lev}.NE))-2*mesh{lev}.NE;
cdE=contour_dof{lev}(find(contour_dof{lev}<2*mesh{lev}.NE));
nodecdN=mesh{lev}.node(cdN,[1,2]);
nodecdE=0.5*(mesh{lev}.node(mesh{lev}.edge(cdE,1),1:2)+mesh{lev}.node(mesh{lev}.edge(cdE,2),1:2));


scatter(nodecdE(:,1),nodecdE(:,2),200,'filled','MarkerFaceColor',[0.7 0.7 0.7]);
scatter(nodecdN(:,1),nodecdN(:,2),200,'filled','MarkerFaceColor',[0. 0. 0.]);

end

end

end