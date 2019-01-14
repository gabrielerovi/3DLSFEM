function print_scalar_function(mesh,f)

L=size(mesh);
L=L(1);
figure
cont_elem=1;
for t=1:mesh{L}.NT
    elem=mesh{L}.elem(t,:);
    node=mesh{L}.node(elem,:);
    for ii=1:length(node(:,1))
        v(ii,1)=f(node(ii,1),node(ii,2));
    end

    fill3(node(:,1),node(:,2),v,v);
    hold on
end







end