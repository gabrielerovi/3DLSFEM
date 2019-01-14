function alpha= BarycentricCoordinatesTriangle(P, nodes_triangleABC)

Area=AreaTriangle(nodes_triangleABC) ;

nodesPBC=nodes_triangleABC;
nodesPBC(1,:)=P;
alpha(1)= AreaTriangle(nodesPBC) /(Area);

nodesPAC=nodes_triangleABC;
nodesPAC(2,:)=P;

alpha(2)=AreaTriangle(nodesPAC) /(Area);

alpha(3)=1-alpha(1)-alpha(2);

for kk=1:3
if(alpha(kk)<0 || alpha(kk)>1)
msg = 'The point p is not inside the triangle';
error(msg)
end
end

end