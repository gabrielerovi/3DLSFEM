
function T_edge_bc=T_to_bc2D(elem,boundary,edge_per_elem)
NT=length(elem(:,1));
NB=length(boundary(:,1));
for t=1:NT
    
    elemT=sort(elem(t,:));
    edge=[elemT(1),elemT(2);elemT(1),elemT(3);elemT(2),elemT(3);] ;
    
    for e=1:edge_per_elem
        
        cont=0;
        for bb=1:NB
            if(sort(boundary(bb,[1,2]))==edge(e,[1,2]) )
                T_edge_bc(t,e)=boundary(bb,3);   
                break;
            end
            cont=cont+1;
        end
                    if(cont==NB)
                     T_edge_bc(t,e)=0;
                    end
    end
end