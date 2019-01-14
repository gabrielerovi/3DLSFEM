function RT_normal =RT03Dnormal(nodes_coord_of_face)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Normal of a face of the mesh whose nodes have coordinates side %%%%%
%%%%% face=mesh.face(ff,:);      (n1<n2<n3, already ordered)         %%%%%
%%%%% nodes_coord_of_face=mesh.node(face,:);                         %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

side=nodes_coord_of_face;
RT_normal=cross(side(2,:)-side(1,:),side(3,:)-side(1,:))';
RT_normal=RT_normal/norm(RT_normal);
end