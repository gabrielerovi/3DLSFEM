function meshes_write(mesh)
L=size(mesh);
L=L(1);
for lev=1:L
    string1='SoniaLaFurby_';
    string2=num2str(lev);
    string3='.txt';
    string=strcat(string1,string2,string3);
    node=mesh{lev}.node;
    elem=mesh{lev}.elem;
    boundary=mesh{lev}.boundary;
    lengths=[length(node(:,1)),length(elem(:,1)),length(boundary(:,1))];
    fileID = fopen(string,'w');
    
    fprintf(fileID,'%d %d %d\n',lengths);
    for nn=1:lengths(1)
    fprintf(fileID,'%f %f\n',node(nn,:));
    end
    for tt=1:lengths(2)
    fprintf(fileID,'%d %d %d\n',elem(tt,:));
    end
    for bb=1:lengths(3)
    fprintf(fileID,'%d %d %d\n',boundary(bb,:));
    end
    fclose(fileID);
end
end