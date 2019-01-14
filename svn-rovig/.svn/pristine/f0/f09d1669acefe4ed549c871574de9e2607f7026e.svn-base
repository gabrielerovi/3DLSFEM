function mesh=mesh_make_it_circle(mesh)

radius=0.5 * (max(mesh{1}.node(:,1))-min(mesh{1}.node(:,1)));
center(1)=0;
center(2)=radius;
L=length(mesh);
for lev=2:L
boundary=mesh{lev}.boundary;
for bb=1:length(boundary(:,1))
for nn=boundary(bb,1:end-1)
    if( mesh{lev}.node(nn,2)<radius)
    val=radius/sqrt( (mesh{lev}.node(nn,1)-center(1))^2 + (mesh{lev}.node(nn,2)-center(2))^2 );
    mesh{lev}.node(nn,1)= val*(mesh{lev}.node(nn,1)-center(1)) + center(1);
    mesh{lev}.node(nn,2)= val*(mesh{lev}.node(nn,2)-center(2)) + center(2);
    end
end
end


end

end