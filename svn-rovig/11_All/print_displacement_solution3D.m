function print_displacement_solution3D(mesh,sol)

L=length(mesh);
NF=mesh{L}.NF;
N=mesh{L}.N;
F_components=mesh{1}.F_components;
figure
increase=1;
cont_elem=1;

disp=sparse(N,3);
disp(:,1)=sol(F_components*NF+1:F_components*NF+N);
disp(:,2)=sol(F_components*NF+N+1:F_components*NF+2*N);
disp(:,3)=sol(F_components*NF+2*N+1:F_components*NF+3*N);

newnode=mesh{L}.node+disp;
% showmesh3(mesh{L}.node,mesh{L}.elem);

hold on

showmesh3(newnode,mesh{L}.elem);


figure
hold on


for ff=1:mesh{L}.NF
    onboundary=mesh{L}.F_bc(ff);
    if(onboundary>0)
    elem=mesh{L}.face(ff,:);
    faces=mesh{L}.node([elem(1),elem(2),elem(3)],:);
    node=mesh{L}.node(elem,:);
    w=[0.2;0.2;0.2];
    %fill3(node(:,1),node(:,2),w,w);
     
    hold on
    h=fill3(faces(:,1),faces(:,2),faces(:,3),'yellow');
    set(h,'facealpha',1);
    loc_disp=full(increase*disp(elem,:));
    
    faces=faces+loc_disp;
    
     
    hold on
    
    h=fill3(faces(:,1),faces(:,2),faces(:,3),'yellow');
    set(h,'facealpha',1);
    
    for ii=1:length(elem)
    quiver3(node(ii,1),node(ii,2),node(ii,3),loc_disp(ii,1),loc_disp(ii,2),loc_disp(ii,3),'-','LineWidth',2,'MaxHeadSize',0.6,'color','red','AutoScale','on', 'AutoScaleFactor', 1)
    end

    end
end







% for t=1:mesh{L}.NT
%     elem=mesh{L}.elem(t,:);
%     faces1=mesh{L}.node([elem(1),elem(2),elem(3)],:);
%     faces2=mesh{L}.node([elem(1),elem(2),elem(4)],:);
%     faces3=mesh{L}.node([elem(2),elem(3),elem(4)],:);
%     faces4=mesh{L}.node([elem(1),elem(3),elem(4)],:);
%     node=mesh{L}.node(elem,:);
%     w=[0.2;0.2;0.2];
%     %fill3(node(:,1),node(:,2),w,w);
%      
%     hold on
%     h=fill3(faces1(:,1),faces1(:,2),faces1(:,3),'red');
%     set(h,'facealpha',.2);
%     h=fill3(faces2(:,1),faces2(:,2),faces2(:,3),'red');
%     set(h,'facealpha',.2);
%     h=fill3(faces3(:,1),faces3(:,2),faces3(:,3),'red');
%     set(h,'facealpha',.2);
%     h=fill3(faces4(:,1),faces4(:,2),faces4(:,3),'red');
%     set(h,'facealpha',.2);
% 
%     
%     loc_disp=full(increase*disp(elem,:));
%     loc_disp1=loc_disp([1,2,3],:);
%     loc_disp2=loc_disp([1,2,4],:);
%     loc_disp3=loc_disp([1,3,4],:);
%     loc_disp4=loc_disp([1,3,4],:);
%     
%     faces1=faces1+loc_disp1;
%     faces2=faces2+loc_disp2;
%     faces3=faces3+loc_disp3;
%     faces4=faces4+loc_disp4;
%     
%     
%     w=[0.5;0.5;0.5];
%     %fill3(node(:,1),node(:,2),w,w);
%      
%     hold on
% %     for ii=1:length(elem)
% %     quiver3(node(ii,1),node(ii,2),node(ii,3),loc_disp(ii,1),loc_disp(ii,2),loc_disp(ii,3),'-','LineWidth',2,'MarkerEdgeColor','r');
% %     end
%     h=fill3(faces1(:,1),faces1(:,2),faces1(:,3),'magenta');
%     set(h,'facealpha',.15);
%     h=fill3(faces2(:,1),faces2(:,2),faces2(:,3),'magenta');
%     set(h,'facealpha',.15);
%     h=fill3(faces3(:,1),faces3(:,2),faces3(:,3),'magenta');
%     set(h,'facealpha',.15);
%     h=fill3(faces4(:,1),faces4(:,2),faces4(:,3),'magenta');
%     set(h,'facealpha',.15);
%     
%     
%     
%     
%     
% end
% 
% 
% 
% % for t=1:mesh{L}.NT
% %     elem=mesh{L}.elem(t,:);
% %     faces1=mesh{L}.node([elem(1),elem(2),elem(3)],:);
% %     faces2=mesh{L}.node([elem(1),elem(2),elem(4)],:);
% %     faces3=mesh{L}.node([elem(2),elem(3),elem(4)],:);
% %     faces4=mesh{L}.node([elem(1),elem(3),elem(4)],:);
% %     node=mesh{L}.node(elem,:);
% %     loc_disp=full(increase*disp(elem,:));
% %     loc_disp1=loc_disp([1,2,3],:);
% %     loc_disp2=loc_disp([1,2,4],:)
% %     loc_disp3=loc_disp([1,3,4],:)
% %     loc_disp4=loc_disp([1,3,4],:)
% %     
% %     faces1=faces1+loc_disp1;
% %     faces2=faces2+loc_disp2;
% %     faces3=faces3+loc_disp3;
% %     faces4=faces4+loc_disp4;
% %     
% %     
% %     w=[1;1;1];
% %     %fill3(node(:,1),node(:,2),w,w);
% %      
% %     hold on
% % %     for ii=1:length(elem)
% % %     quiver3(node(ii,1),node(ii,2),node(ii,3),loc_disp(ii,1),loc_disp(ii,2),loc_disp(ii,3),'-','LineWidth',2,'MarkerEdgeColor','r');
% % %     end
% %     h=fill3(faces1(:,1),faces1(:,2),faces1(:,3),w);
% %     set(h,'facealpha',.5);
% %     h=fill3(faces2(:,1),faces2(:,2),faces2(:,3),w);
% %     set(h,'facealpha',.5);
% %     h=fill3(faces3(:,1),faces3(:,2),faces3(:,3),w);
% %     set(h,'facealpha',.5);
% %     h=fill3(faces4(:,1),faces4(:,2),faces4(:,3),w);
% %     set(h,'facealpha',.5);
% % 
% % end
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % % for t=1:mesh{L}.NT
% % %     elem=mesh{L}.elem(t,:);
% % %     faces1=mesh{L}.node([elem(1),elem(2),elem(3)],:);
% % %     faces2=mesh{L}.node([elem(1),elem(2),elem(4)],:);
% % %     faces3=mesh{L}.node([elem(2),elem(3),elem(4)],:);
% %     faces4=mesh{L}.node([elem(1),elem(3),elem(4)],:);
% %     node=mesh{L}.node(elem,:);
% %     loc_disp=increase*disp(elem,:);
% %     w=[0.1;0.1;0.1];
% %     %fill3(node(:,1),node(:,2),w,w);
% %      
% %     hold on
% %     for ii=1:length(elem)
% %     quiver3(node(ii,1),node(ii,2),node(ii,3),loc_disp(ii,1),loc_disp(ii,2),loc_disp(ii,3),'-','LineWidth',2,'MarkerEdgeColor','r');
% %     end
% %     h=fill3(faces1(:,1),faces1(:,2),faces1(:,3),w);
% %     set(h,'facealpha',.5);
% %     h=fill3(faces2(:,1),faces2(:,2),faces2(:,3),w);
% %     set(h,'facealpha',.5);
% %     h=fill3(faces3(:,1),faces3(:,2),faces3(:,3),w);
% %     set(h,'facealpha',.5);
% %     h=fill3(faces4(:,1),faces4(:,2),faces4(:,3),w);
% %     set(h,'facealpha',.5);
% % 
% % end
% % 
% 
% 
% % for t=1:mesh{L}.NT
% %     elem=mesh{L}.elem(t,:);
% %     node=mesh{L}.node(elem,:);
% %     loc_disp1=d1(elem)';
% %     loc_disp2=d2(elem)';
% %     w=[0;0;0];
% %     %fill3(node(:,1),node(:,2),w,w);
% %     new_node(:,1)=node(:,1)+increase*loc_disp1;
% %     new_node(:,2)=node(:,2)+increase*loc_disp2;
% %     
% %     h=fill3(new_node(:,1),new_node(:,2),w,w);
% %     set(h,'facealpha',.5)
%    quiver(new_node(:,1),new_node(:,2),loc_disp1,loc_disp2,'-','LineWidth',2,'MarkerEdgeColor','r');
% end
% 





end