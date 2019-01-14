function print_scalar_solution(mesh,x)

L=size(mesh);
L=L(1);
figure
h=colorbar('northoutside');
t=get(h,'Limits');
set(h,'Ticks',linspace(t(1),t(2),5));
set(gca,'FontSize',20);
colormap('jet');
cont_elem=1;
minimum=0;
maximum=0;
for t=1:mesh{L}.NT
    elem=mesh{L}.elem(t,:);
    node=mesh{L}.node(elem,:);
    v=x(elem);
    ss=size(v);
    if(ss(1)<ss(2))
        v=v';
    end
    fill3(node(:,1),node(:,2),v,v);
    minimum=(min(minimum,min(v)));
    maximum=(max(maximum,max(v)));
    hold on
end

caxis([minimum,maximum]);





end