function [NF,face,elemF,F_to_T]=create_face_structures3D(face_per_elem,node,elem)
% face_per_elem =3/4 (2D/3D for triangles/tetrahedra)
% in 2D we do not use faces, but directly edges

% N = number of vertices
% NT = number of elements
% NF = the number of faces

% node = matrix with nodes coordinates
% elem = matrix with nodes dofs for each element (raw)
% elemF = matrix with face dofs for each element (raw)


N=length(node(:,1));
NT=length(elem(:,1));
NF=0;
%create a sufficient big cell
F_to_T=cell(NT*4,1);
%face connectivity matrix. 1 face exists if the two nodes n1 and n2 are
%connected
face_connectivity=zeros(N,N,N);

for T=1:NT
    
   vertex=sort(elem(T,:));
   faces(1,:)=[vertex(1),vertex(2),vertex(3)];
   faces(2,:)=[vertex(1),vertex(2),vertex(4)];
   faces(3,:)=[vertex(1),vertex(3),vertex(4)];
   faces(4,:)=[vertex(2),vertex(3),vertex(4)];
   
   
   % in 2D/3D, there are 3/4 faces (but we use this function only for 3D)
   for n=1:face_per_elem
       % vertices n1 and n2, with n2 >n1
       n1=faces(n,1);
       n2=faces(n,2);
       n3=faces(n,3);
       
   if(face_connectivity( n1, n2,  n3 )==0)  
       % +1 edge dof
       NF=NF+1;
       % update edge_connectivity
       face_connectivity( n1, n2, n3 ) = NF;
       face_connectivity( n1, n2, n3 ) = NF;
       % update edge with the node dofs 
       face(NF,[1,2,3])=faces(n,[1,2,3]);
   end
   
   NF_tmp=face_connectivity( n1, n2, n3 );
   
   if(isempty(F_to_T{NF_tmp}))
   F_to_T{NF_tmp}=cell(1,1);
   F_to_T{NF_tmp}{end}=T;
   else
   F_to_T{NF_tmp}{end+1}=T;   
   end
   % create a structure like elem, but with edge dofs instead of vertex dofs 
   elemF(T,n)=NF_tmp; 
   end

end

 F_to_T=F_to_T(~cellfun('isempty',F_to_T))  ;
end

