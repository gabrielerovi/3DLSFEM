function print_vector_solution(mesh,x)

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


for t=1:mesh{L}.NT
    elem=mesh{L}.elem(t,:);
    elemE=mesh{L}.elemE(t,:);
    node=mesh{L}.node(elem,:);
    v=x(elemE);
    phi=phiRT2D(node,node);
    vector=zeros(length(phi(:,1)),2);
    
    barycenter=sum(node)/length(node(:,1));
    for e=1:length(phi(:,1))
        for i=1:3
            for dd=1:2
            vector(e,dd)=vector(e,dd)+v(i)*phi(i,e,dd);
            end
        end
        %quiver(node(e,1),node(e,2),vector(e,1),vector(e,2),'-','LineWidth',2,'MarkerEdgeColor','r');
    end
    quiver(barycenter(1),barycenter(2),vector(e,1),vector(e,2),'-','LineWidth',2,'MarkerEdgeColor','r');

end

end