function [node,elem,boundary]=mesh_read3D(input_file)

    file_id=fopen(input_file);
    for ii=1:5
    data_line=fgetl(file_id);
    end 
    N=sscanf(data_line,'%i');
    node=[];
    for ii=1:N
        data_line=fgetl(file_id);
        node_coord = sscanf(data_line,'%f')';
        node=[node; node_coord(2:end)];
    end

     for ii=1:3
    data_line=fgetl(file_id);
    end 

    NTot=sscanf(data_line,'%i');
    elem=[];
    boundary=[];
    
    cont=0;
    old_boundary=0;
     for ii=1:NTot
        data_line=fgetl(file_id);
        elem_dof = sscanf(data_line,'%f')';
        %% it 
        if(length(elem_dof)<8)
        % then it is an edge
        elseif(length(elem_dof)<9)
        % then it is boundary elem
        if(old_boundary~=elem_dof(4))
            old_boundary=elem_dof(4);
            cont=cont+1;
        end
        boundary=[boundary; elem_dof(end-2:end),cont];
        elseif(length(elem_dof)<10)
        % then it is boundary elem
        elem=[elem; elem_dof(end-3:end)];
        end
     end
    
     
   fclose(file_id)

end