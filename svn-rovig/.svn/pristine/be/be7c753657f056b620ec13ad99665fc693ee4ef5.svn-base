function print_stress_solution(mesh,x,component)

L=size(mesh);
L=L(1);
figure
hold on
for t=1:mesh{L}.NT
    elem=mesh{L}.elem(t,:);
    elemE=mesh{L}.elemE(t,:);
    node=mesh{L}.node(elem,:);
    line(node([1:end,1],1),node([1:end,1],2));
end

maximum=0;
for t=1:mesh{L}.NT
    elem=mesh{L}.elem(t,:);
    elemE=mesh{L}.elemE(t,:);
    node=mesh{L}.node(elem,:);
    v=x(elemE);
    phi=phiRT2D(node,node);
    vector=zeros(length(phi(:,1)),2);
    nodeE= [ 0.5 * ( node(1,1)+node(2,1) ),0.5 * (node(1,2)+node(2,2));
             0.5 * ( node(2,1)+node(3,1) ),0.5 * (node(2,2)+node(3,2));
             0.5 * ( node(1,1)+node(3,1) ),0.5 * (node(1,2)+node(3,2));];
    
    barycenter=sum(node)/length(node(:,1));
    for e=1:length(phi(:,1))
        for i=1:3
            for dd=1:2
            vector(e,dd)=vector(e,dd)+v(i)*phi(i,e,dd);
            end
        end
       % scatter3(nodeE(e,1),nodeE(e,2),vector(e,component),100.0,[0,0,0]);
      %  quiver(node(e,1),node(e,2),vector(e,1),vector(e,2),'-','LineWidth',2,'MarkerEdgeColor','r');
    end
   %maximum=max(maximum,max(abs(vector(:,component))));

     v=[vector(e,component);vector(e,component);vector(e,component)];
     fill3(node(:,1),node(:,2),v,v);
% 
end
% 
% for t=1:mesh{L}.NT
%     elem=mesh{L}.elem(t,:);
%     elemE=mesh{L}.elemE(t,:);
%     node=mesh{L}.node(elem,:);
%     v=x(elemE);
%     vector=zeros(length(phi(:,1)),2);
%     nodeE= [ 0.5 * ( node(1,1)+node(2,1) ),0.5 * (node(1,2)+node(2,2));
%              0.5 * ( node(2,1)+node(3,1) ),0.5 * (node(2,2)+node(3,2));
%              0.5 * ( node(1,1)+node(3,1) ),0.5 * (node(1,2)+node(3,2));];
%     phi=phiRT2D(nodeE,node);
%     
%     barycenter=sum(node)/length(node(:,1));
%     for e=1:length(phi(:,1))
%             %for dd=1:2
%             %vector(e,dd)=vector(e,dd)+v(i)*phi(i,e,dd);
%             vector(e)=v(e)* (phi(e,e,1) * 0 +phi(e,e,2) * 1);
%         
%         val=abs(vector(e,component))/maximum;
%         if(abs(vector(e,component))==0.0 || abs(vector(e,component)-(-0.0005))<0.0001)
%         scatter3(nodeE(e,1),nodeE(e,2),vector(e),100.0,[val,val,val],'filled');
%         end
%         %quiver(node(e,1),node(e,2),vector(e,1),vector(e,2),'-','LineWidth',2,'MarkerEdgeColor','r');
%     end
% 
% %     v=[vector(e,component);vector(e,component);vector(e,component)];
% %     fill3(node(:,1),node(:,2),v,v);
% % 
% end
% 


end