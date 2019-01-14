function [node,elem,boundary]=mesh_read(input_file)

    file_id=fopen(input_file);
    data_line=fgetl(file_id);
    data_mesh=sscanf(data_line,'%i');
    N=data_mesh(1);
    NT=data_mesh(2);
    NB=data_mesh(3);
    
    node=[];
    for ii=1:N
        data_line=fgetl(file_id);
        node_coord = sscanf(data_line,'%f')';
        node=[node; node_coord];
    end
    
    elem=[];
     for ii=1:NT
        data_line=fgetl(file_id);
        elem_dof = sscanf(data_line,'%f')';
        elem=[elem; elem_dof];
     end
    
    boundary=[];
     for ii=1:NB
        data_line=fgetl(file_id);
        elem_boundary = sscanf(data_line,'%f')';
        boundary=[boundary; elem_boundary];
     end
     
     
   fclose(file_id)

end