function phidotn=phi_dot_n(mesh,ff)



face=mesh.face(ff,:);
side=mesh.node(face,:);
normal=mesh.normal_face{ff};
RT_normal =RT03Dnormal(side);
%RT_normal=cross(side(2,:)-side(1,:),side(3,:)-side(1,:))';

if(RT_normal'*normal>0)
    sign=1;
else
    sign=-1;
end

A=AreaTriangle(side);

phidotn=sign/A;
        
end