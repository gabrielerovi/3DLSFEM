function [bool,Marker]=is_surface_dirichlet(label,type_of_dof)



[dirichlet,n_and_or_t,bool_bc]= boundary_value_bool(type_of_dof);



cont=0;
Marker=[];
% dof=node
if(type_of_dof==1)
    for i=1:length(label)
        if(label(i)>0) 
        bool(i) = bool_bc(label(i));
        if(bool(i)>0)
        cont=cont+1;
        Marker(cont)=label(i);
        end
        else
        bool(i) = 0;
        end
        
    end
elseif(type_of_dof==2)
    for i=1:length(label)
        if(label(i)>0)   
                    
        bool(i) = bool_bc(label(i));
        if(bool(i)>0)
        cont=cont+1;
        Marker(cont)=label(i); 
        end
        else
        bool(i) = 0;
        end
    end
elseif(type_of_dof==3)
    for i=1:length(label)
        if(label(i)>0)             
        bool(i) = bool_bc(label(i));
        if(bool(i)>0)
        cont=cont+1; 
        Marker(cont)=label(i);
        end
        else
        bool(i) = 0;
        end
    end
end

end