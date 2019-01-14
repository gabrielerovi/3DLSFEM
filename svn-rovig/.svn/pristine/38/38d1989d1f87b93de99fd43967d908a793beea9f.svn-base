function meshwidth=hminhmax(mesh)
hmin=1000000;
hmax=0;

% in 2D
if(mesh.face_per_elem==3)
dim=2;
NF=mesh.NE;
face=mesh.edge;
else
dim=3;
NF=mesh.NF;
face=mesh.face;
end
node=mesh.node;

for ff=1:NF
    
    nodeff=node(face(ff,:),1:dim);
    
    for nn=1:length(nodeff(:,1))-1
        h=norm(nodeff(nn,:)-nodeff(nn+1,:));
        hmin=min(hmin,h);
        hmax=max(hmax,h);
    end
        h=norm(nodeff(1,:)-nodeff(end,:));
        hmin=min(hmin,h);
        hmax=max(hmax,h);
      
        
end

meshwidth.hmin=hmin;
meshwidth.hmax=hmax;
meshwidth.hmaxfrachmin=hmax/hmin;
end