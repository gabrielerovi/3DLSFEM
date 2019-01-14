function phidotn=phi_dot_n(mesh,ee)



edge=mesh.edge(ee,:);
side=mesh.node(edge,:);
normal=mesh.normal_edge{ee};
RT_normal=[ side(2,2)-side(1,2);
           -side(2,1)+side(1,1); ];
if(RT_normal'*normal>0)
    sign=1;
else
    sign=-1;
end
phidotn=sign/norm(side(2,[1,2])-side(1,[1,2]));
        
end