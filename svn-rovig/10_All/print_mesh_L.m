function print_mesh_L(mesh,uwant2contour,contour_dof)
L=length(mesh);
for lev=1:L
figure
for ii=1:mesh{lev}.N
    
x = mesh{lev}.node(ii,1);
y = mesh{lev}.node(ii,2);
z = mesh{lev}.node(ii,3);
txt = text(x,y,z,num2str(ii));
txt.Color='r';  
end

for ii=1:mesh{lev}.NF
side=mesh{lev}.face(ii,:);
node1=mesh{lev}.node(side(1),:);
node2=mesh{lev}.node(side(2),:);
node3=mesh{lev}.node(side(3),:);


coeff1=1;%NDtoRT{L}(ii,mesh{L}.edge(ii,1));
coeff2=1;%NDtoRT{L}(ii,mesh{L}.edge(ii,2));

x = 1/3 * (node1(1)+node2(1)+node3(1));
y = 1/3 * (node1(2)+node2(2)+node3(2));
z = 1/3 * (node1(3)+node2(3)+node3(3));
txt = text(x,y,z,num2str(ii)) ;
txt.Color='b';  

hold on
% Their vertial concatenation is what you want
pts1 = [node1(1),node1(2),node1(3); node2(1),node2(2),node2(3)];
pts2 = [node1(1),node1(2),node1(3); node3(1),node3(2),node3(3)];
pts3 = [node2(1),node2(2),node2(3); node3(1),node3(2),node3(3)];

% Because that's what line() wants to see    
line(pts1(:,1), pts1(:,2), pts1(:,3))
line(pts2(:,1), pts2(:,2), pts2(:,3))
line(pts3(:,1), pts3(:,2), pts3(:,3))
% % Alternatively, you could use plot3:
% plot3(pts(:,1), pts(:,2), pts(:,3))
% 
% %quiver(node1(1),node1(2),node1(3),arrow(1),arrow(2),arrow(3),0)
% plot(plot::Line3d([node1(1),node1(2),node1(3)],[node2(1),node2(2),node2(3)]));
% plot(plot::Line3d([node1(1),node1(2),node1(3)],[node3(1),node3(2),node3(3)]));
% plot(plot::Line3d([node2(1),node2(2),node2(3)],[node2(1),node2(2),node2(3)]));
end

if(uwant2contour{lev})
cdN=contour_dof{lev}(find(contour_dof{lev}>2*mesh{lev}.NF))-2*mesh{lev}.NF;
cdF=contour_dof{lev}(find(contour_dof{lev}<2*mesh{lev}.NF));
nodecdN=mesh{lev}.node(cdN,[1,2]);
nodecdF=0.5*(mesh{lev}.node(mesh{lev}.face(cdF,1),1:2)+mesh{lev}.node(mesh{lev}.face(cdF,2),1:2));


scatter(nodecdF(:,1),nodecdE(:,2),200,'filled','MarkerFaceColor',[0.7 0.7 0.7]);
scatter(nodecdN(:,1),nodecdN(:,2),200,'filled','MarkerFaceColor',[0. 0. 0.]);

end

end

end