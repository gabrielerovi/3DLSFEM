
function T_face_bc=T_to_bc3D(elem,boundary,face_per_elem)
NT=length(elem(:,1));
NB=length(boundary(:,1));
for t=1:NT
    
    elemT=sort(elem(t,:));
    face=[elemT(1),elemT(2),elemT(3);
          elemT(1),elemT(2),elemT(4);
          elemT(1),elemT(3),elemT(4);
          elemT(2),elemT(3),elemT(4);] ;
    
    for ff=1:face_per_elem
        
        cont=0;
        for bb=1:NB
            if(sort(boundary(bb,[1,2,3]))==face(ff,[1,2,3]) )
                T_face_bc(t,ff)=boundary(bb,4);   
                break;
            end
            cont=cont+1;
        end
                    if(cont==NB)
                     T_face_bc(t,ff)=0;
                    end
    end
end