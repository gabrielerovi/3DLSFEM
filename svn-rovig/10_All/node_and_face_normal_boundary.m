function [mesh]=node_and_face_normal_boundary(mesh,parameters)
L=length(mesh);
kk=[3;1;2];
for lev=1:L
grid=mesh{lev};
node_per_elem=grid.node_per_elem;
face_per_elem=grid.face_per_elem;
node=grid.node;
face=grid.face;
NT=grid.NT;
N=grid.N;
NF=grid.NF;
F_bc=grid.F_bc;

[dirichlet_F,n_and_or_t_F,bool_bc_F]= boundary_value_bool(3);



for ff=1:NF
    normal_face{lev,ff}=[];
end

for tt=1:NT
elemF=grid.elemF(tt,:);
elem=grid.elem(tt,:);


for ff=1:face_per_elem
    fftot=elemF(ff);
%     if(F_bc(fftot)>0)
        vertices=face(fftot,:);
        opposite_vertex=setdiff(elem,face(fftot,:));
        side=node(vertices,:);
        opposite_node=node(opposite_vertex,:);
        mean_node=mean(side);
        opposite_vector=mean_node-opposite_node;
        normal_face{lev,fftot}=cross(side(2,:)-side(1,:),side(3,:)-side(1,:))';
        normal_face{lev,fftot}=sign((opposite_vector*normal_face{lev,fftot}))*normal_face{lev,fftot}/norm(normal_face{lev,fftot});       
%     end
end
end

end


for lev=1:L
    for ff=1:mesh{lev}.NF
    mesh{lev}.normal_face{ff}=normal_face{lev,ff};
    end
end
end