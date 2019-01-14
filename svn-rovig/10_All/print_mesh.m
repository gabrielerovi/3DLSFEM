function print_mesh(mesh,dim)
figure
hold on
L=size(mesh);
L=L(1);
if(dim==2)
 for nn=1:mesh{L}.N
    scatter(mesh{L}.node(nn,1),mesh{L}.node(nn,2),'k','fill');
end
for bb=1:length( mesh{L}.boundary(:,1))
    n1=mesh{L}.node( mesh{L}.boundary(bb,1),:);
    n2=mesh{L}.node( mesh{L}.boundary(bb,2),:);
    
    quiver(n1(1),n1(2), n2(1)-n1(1),n2(2)-n1(2),'filled','r');
end
else
    for nn=1:mesh{L}.N
    scatter3(mesh{L}.node(nn,1),mesh{L}.node(nn,2),mesh{L}.node(nn,3),'k','filled');
    end
  for bb=1:length( mesh{L}.boundary(:,1))
    n1=mesh{L}.node( mesh{L}.boundary(bb,1),:);
    n2=mesh{L}.node( mesh{L}.boundary(bb,2),:);
    n3=mesh{L}.node( mesh{L}.boundary(bb,3),:);
    X=[n1(1);n2(1);n3(1)];
    Y=[n1(2);n2(2);n3(2)];
    Z=[n1(3);n2(3);n3(3)];
    C=[mesh{L}.boundary(bb,4)/50.0 mesh{L}.boundary(bb,4)/50.0 mesh{L}.boundary(bb,4)/50.0];
    hold on
    fill3(X,Y,Z,C);
  end 
figure
showmesh3(mesh{L}.node,mesh{L}.elem)
end
if(dim==2)
 for nn=1:mesh{L}.N
    scatter(mesh{L}.node(nn,1),mesh{L}.node(nn,2),'k','fill');
end
for bb=1:length( mesh{L}.boundary(:,1))
    n1=mesh{L}.node( mesh{L}.boundary(bb,1),:);
    n2=mesh{L}.node( mesh{L}.boundary(bb,2),:);
    
    quiver(n1(1),n1(2), n2(1)-n1(1),n2(2)-n1(2),'filled','r');
end
else
    for nn=1:mesh{L}.N
    scatter3(mesh{L}.node(nn,1),mesh{L}.node(nn,2),mesh{L}.node(nn,3),'k','filled');
    end
  for bb=1:length( mesh{L}.boundary(:,1))
    n1=mesh{L}.node( mesh{L}.boundary(bb,1),:);
    n2=mesh{L}.node( mesh{L}.boundary(bb,2),:);
    n3=mesh{L}.node( mesh{L}.boundary(bb,3),:);
    X=[n1(1);n2(1);n3(1)];
    Y=[n1(2);n2(2);n3(2)];
    Z=[n1(3);n2(3);n3(3)];
    C=[mesh{L}.boundary(bb,4)/50.0 mesh{L}.boundary(bb,4)/50.0 mesh{L}.boundary(bb,4)/50.0];
    hold on
    fill3(X,Y,Z,C);
  end 
figure
showmesh3(mesh{L}.node,mesh{L}.elem);
end

end