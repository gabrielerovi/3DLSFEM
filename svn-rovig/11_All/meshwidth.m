function h=meshwidth(mesh)
L=size(mesh);
L=L(1);

for lev=1:L
    grid=mesh{lev};
    NE=mesh{lev}.NE;
    h(lev)=10000000;
    for ee=1:NE
        vertices=grid.edge(ee,:);
        node=grid.node(vertices,:);
        h(lev)=min(h(lev),norm(node(1,:)-node(2,:)));
    end
    
end



end