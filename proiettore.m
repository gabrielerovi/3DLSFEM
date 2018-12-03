close all
clear all
clc

L=2;
coordinates_node_='coordinates_node_'
coordinates_edge_='coordinates_edge_'
dofmap_edge_='dofmap_edge_'
dofmap_node_='dofmap_node_'
load('solutioncoarse.mat')
load('solutionfine.mat')
%load('coordinates_coarse.mat')
load('coordinates_edge_coarse.mat')
load('coordinates_edge_fine.mat')
load('coordinates_fine.mat')

load('dofmap_edge_coarse.mat')
load('dofmap_edge_fine.mat')
load('dofmap_node_coarse.mat')
load('dofmap_node_fine.mat')
load('WC2WFnumpy.mat')
% we are in matlab, not c++
dofmap_edge_coarse=dofmap_edge_coarse+1;
dofmap_edge_fine=dofmap_edge_fine+1;
dofmap_node_coarse=dofmap_node_coarse+1;
dofmap_node_fine=dofmap_node_fine+1;



for lev=0:L-1
    load(strcat(coordinates_node_,num2str(lev),'.mat'))
    coordinates_node_list{lev+1}=coordinates_node;
    load(strcat(coordinates_edge_,num2str(lev),'.mat'))
    coordinates_edge_list{lev+1}=coordinates_edge;
    load(strcat(dofmap_node_,num2str(lev),'.mat'))
    dofmap_node_list{lev+1}=dofmap_node;
    load(strcat(dofmap_edge_,num2str(lev),'.mat'))
    dofmap_edge_list{lev+1}=dofmap_edge;
    
    
    figure
    hold on

    for nn=1:length(coordinates_node_list{lev+1})
    x=coordinates_node_list{lev+1}(nn,1);
    y=coordinates_node_list{lev+1}(nn,2);
    z=coordinates_node_list{lev+1}(nn,3);
    %solc_loc=[solc(dofmap_node_coarse(nn,1)),solc(dofmap_node_coarse(nn,2))];
    %scatter(x,y,'k','fill'); 
    str=strcat(num2str(nn-1),':__',num2str(dofmap_node_list{lev+1}(nn,:)))%,':__',num2str(solc_loc));
    txt = text(x,y,z,str) ;
    txt.Color='b';
    end
    
    for ee=1:length(coordinates_edge_list{lev+1})
    x=coordinates_edge_list{lev+1}(ee,1);
    y=coordinates_edge_list{lev+1}(ee,2);
    z=coordinates_edge_list{lev+1}(ee,3);
    %solc_loc=[solc(dofmap_edge_coarse(ee,1)),solc(dofmap_edge_coarse(ee,2))];
    
    str=strcat(num2str(ee-1),':__',num2str(dofmap_edge_list{lev+1}(ee,:)))%,':__',num2str(solc_loc));
    txt = text(x,y,z,str) ;
    txt.Color='r';
    end

    
end






for nn=1:length(coordinates_coarse)
    x=coordinates_coarse(nn,1);
    y=coordinates_coarse(nn,2);
    solc_loc=[solc(dofmap_node_coarse(nn,1)),solc(dofmap_node_coarse(nn,2))];
    %scatter(x,y,'k','fill'); 
    str=strcat(num2str(nn),':__',num2str(dofmap_node_coarse(nn,:)),':__',num2str(solc_loc));
    txt = text(x,y,str) ;
    txt.Color='b';
end

for ee=1:length(coordinates_edge_coarse)
    x=coordinates_edge_coarse(ee,1);
    y=coordinates_edge_coarse(ee,2);
    solc_loc=[solc(dofmap_edge_coarse(ee,1)),solc(dofmap_edge_coarse(ee,2))];
    
    str=strcat(num2str(ee),':__',num2str(dofmap_edge_coarse(ee,:)),':__',num2str(solc_loc));
    txt = text(x,y,str) ;
    txt.Color='r';
end


figure

for nn=1:length(coordinates_fine)
    x=coordinates_fine(nn,1);
    y=coordinates_fine(nn,2);
    solf_loc=[solf(dofmap_node_fine(nn,1)),solf(dofmap_node_fine(nn,2))];
    str=strcat(num2str(nn),':__',num2str(dofmap_node_fine(nn,:)),':__',num2str(solf_loc));
    txt = text(x,y,str) ;
    txt.Color='b';
end

for ee=1:length(coordinates_edge_fine)
    x=coordinates_edge_fine(ee,1);
    y=coordinates_edge_fine(ee,2);
    solf_loc=[solf(dofmap_edge_fine(ee,1)),solf(dofmap_edge_fine(ee,2))];
    str=strcat(num2str(ee),':__',num2str(dofmap_edge_fine(ee,:)),':__',num2str(solf_loc));
    txt = text(x,y,str) ;
    txt.Color='r';
end




 for nn=1:mesh{L}.N
    scatter(mesh{L}.node(nn,1),mesh{L}.node(nn,2),'k','fill');
end
for bb=1:length( mesh{L}.boundary(:,1))
    n1=mesh{L}.node( mesh{L}.boundary(bb,1),:);
    n2=mesh{L}.node( mesh{L}.boundary(bb,2),:);
    
    %quiver(n1(1),n1(2), n2(1)-n1(1),n2(2)-n1(2),'filled','r');
end

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

