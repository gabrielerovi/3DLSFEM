function savevtkdisplacement(mesh,sol)

fid = fopen('displacement.vtk','wt');
fprintf(fid, '# vtk DataFile Version 3.0\n');
fprintf(fid, 'vtk output\n');
fprintf(fid, 'ASCII\n');
fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
NF=mesh.NF;
N=mesh.N;

disp_x=sol(NF*3+1:NF*3+N);
disp_y=sol(NF*3+N+1:NF*3+2*N);
disp_z=sol(NF*3+2*N+1:NF*3+3*N);
NT=mesh.NT;
node=mesh.node;
elem=mesh.elem;
newline=N;
formatSpec = 'POINTS %d float\n';
fprintf(fid,formatSpec,newline)
for ii=1:N
formatSpec ='%g %g %g\n';
fprintf(fid,formatSpec,node(ii,:))    
end
newline=[NT NT*(4+1)];
formatSpec = 'CELLS %d %d\n';
fprintf(fid,formatSpec,newline);
for ii=1:NT
formatSpec ='%d %d %d %d %d\n';
locelem=elem(ii,:)-1;
fprintf(fid,formatSpec,[4, locelem])    
end

newline=NT;
formatSpec = 'CELL_TYPES %d\n';
fprintf(fid,formatSpec,newline);
for ii=1:NT
fprintf(fid, '10\n');
end

newline=N;
formatSpec = 'POINT_DATA %d\n';
fprintf(fid,formatSpec,newline);

newline=N;
formatSpec = 'VECTORS disp float\n';
fprintf(fid,formatSpec);
for ii=1:N
newline=[disp_x(ii),disp_y(ii),disp_z(ii)];
formatSpec ='%g %g %g\n';
fprintf(fid,formatSpec,newline);  
end

fclose(fid);
end