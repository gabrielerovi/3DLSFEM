function max_err=print_displacement_error(mesh,d,exact)
L=size(mesh);
L=L(1);
figure
increase=1;
cont_elem=1;
max_err=0;
hold on
for t=1:mesh{L}.NT
    elem=mesh{L}.elem(t,:);
    node=mesh{L}.node(elem,:);
    loc_disp=d(elem)';
    for i=1:3
    ex(i,1)=exact(node(i,1),node(i,2));
    end
    err=ex-loc_disp;
    %err=ex;
    fill3(node(:,1),node(:,2),err,err);
    max_err=max(max_err,max(err));
end





set(gca,'FontSize',20)


end