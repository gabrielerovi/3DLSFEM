function [flag_for_bc_face,flag_for_each_face] =flag_boundary_faces3D  (boundary,face)

NF=length(face(:,1));
BF=length(boundary(:,1));

flag_for_bc_face=zeros(BF,1);
flag_for_each_face=zeros(NF,1);
cont=0;

for ff=1:NF
    
    for bb=1:BF
        
        if(sort(boundary(bb,[1,2,3])) == face (ff,:))
            cont=cont+1;
            flag_for_bc_face(cont)=boundary(bb,4)
            flag_for_each_face(ff)=boundary(bb,4);
        end
    end
    
end

end