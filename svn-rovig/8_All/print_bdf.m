function printa_bdf(mesh)

L=size(mesh);
L=L(1);

for lev=1:L
N=mesh{lev}.N;
NT=mesh{lev}.NT;
NB=length(mesh{lev}.boundary(:,1));

node=mesh{lev}.node;
elem=mesh{lev}.elem;
boundary=mesh{lev}.boundary;

format short e
for nn=1:N
    
    Nodes{nn}=['GRID    ',num2str(nn)];
    for ii=length(Nodes{nn})+1:17
         Nodes{nn}(ii)=' ';
    end
    
    Nodes{nn}(17)=num2str(0);
    Nodes{nn}=[ Nodes{nn},'       ', num2str(node(nn,1),'%6.2E\t\t'),num2str(node(nn,2),'%6.2E\t\t'), num2str(0,'%6.2E\t\t')];
        
end


startbb=N+10;
for bb=1:NB
    
    Boundaries{bb}=['CBAR    ',num2str(bb+startbb),...
                    '      ', num2str(boundary(bb,3)),...
                    '       ', num2str(boundary(bb,1)),...
                    '       ', num2str(boundary(bb,2)),...
                    '       ', '0.',...
                    '      ', '0.',...
                    '      ', '0.'];
end

startee=startbb+NB ; %ceil((startbb+NB)/10)*10;
for ee=1:NT
    Elements{ee}=['CTRIA3  ' , num2str(ee+startee),...
                  '      ' , num2str(6),...
                  '       ' , num2str(elem(ee,1)),...
                  '       ' , num2str(elem(ee,2)),...
                  '       ' , num2str(elem(ee,3))];
end




    string1='MeshBDF_';
    string2=num2str(lev);
    string3='.bdf';
    string=strcat(string1,string2,string3);
    fileID = fopen(string,'w');
    
    fprintf(fileID,'%s','$ Created by Gmsh');
    fprintf(fileID,'\n');
    fprintf(fileID,'%s','BEGIN BULK');
    fprintf(fileID,'\n');
    
    for nn=1:N
    fprintf(fileID,'%s %s %s1.2 %s1.2 %s1.2',Nodes{nn});
    fprintf(fileID,'\n');
    end
    
    for bb=1:NB
    fprintf(fileID,'%s %s %s %s %s %s %s',Boundaries{bb});
    fprintf(fileID,'\n');
    end
    
    
    for tt=1:NT
            fprintf(fileID,'%s %s %s %s %s',Elements{tt});
    fprintf(fileID,'\n');
    end
    
    fprintf(fileID,'%s','ENDDATA');
    
    fclose(fileID);

end

end