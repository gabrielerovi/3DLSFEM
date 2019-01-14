function print_displacement_solution(mesh,d1,d2)

L=size(mesh);
L=L(1);
figure
increase=1;
cont_elem=1;

hold on
for t=1:mesh{L}.NT
    elem=mesh{L}.elem(t,:);
    node=mesh{L}.node(elem,:);
    loc_disp1=d1(elem)';
    loc_disp2=d2(elem)';
    w=[0;0;0];
    %fill3(node(:,1),node(:,2),w,w);
    new_node(:,1)=node(:,1);
    new_node(:,2)=node(:,2);
    
    h=fill3(new_node(:,1),new_node(:,2),w,w);
    set(h,'facealpha',.5)
   quiver(new_node(:,1),new_node(:,2),loc_disp1,loc_disp2,'-','LineWidth',2,'MarkerEdgeColor','r');
end



for t=1:mesh{L}.NT
    elem=mesh{L}.elem(t,:);
    node=mesh{L}.node(elem,:);
    loc_disp1=d1(elem)';
    loc_disp2=d2(elem)';
    w=[0;0;0];
    %fill3(node(:,1),node(:,2),w,w);
    new_node(:,1)=node(:,1)+increase*loc_disp1;
    new_node(:,2)=node(:,2)+increase*loc_disp2;
    
    h=fill3(new_node(:,1),new_node(:,2),w,w);
    set(h,'facealpha',.5)
   quiver(new_node(:,1),new_node(:,2),loc_disp1,loc_disp2,'-','LineWidth',2,'MarkerEdgeColor','r');
end






end