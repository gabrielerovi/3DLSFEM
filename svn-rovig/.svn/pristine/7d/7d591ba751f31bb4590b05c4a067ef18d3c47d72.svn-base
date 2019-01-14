function [NE,edge,elemE,E_to_T,N_to_T]=create_edge_structures2D(edge_per_elem,node,elem)
% edge_per_elem =3/6 (2D/3D for triangles/tetrahedra)

% N = number of vertices
% NT = number of elements
% NE = the number of edge

% node = matrix with nodes coordinates
% elem = matrix with nodes dofs for each element (raw)
% elemE = matrix with edge dofs for each element (raw)


N=length(node(:,1));
NT=length(elem(:,1));
NE=0;
node_per_elem=length(elem(1,:));
%create a sufficient big cell
E_to_T=cell(N*N,1);
N_to_T=cell(N,1);
%edge connectivity matrix. 1 edge exists if the two nodes n1 and n2 are
%connected
edge_connectivity=zeros(N,N);

for T=1:NT
    
   vertex=sort(elem(T,:));
   edges=sort([vertex(1,1:end); [vertex(1,2:end),vertex(1,1)]  ])';
   
   % in 2D/3D, there are 3/6 edges
   for n=1:edge_per_elem
       
       % vertices n1 and n2, with n2 >n1
       n1=edges(n,1);
       n2=edges(n,2);
   if(edge_connectivity( n1, n2 )==0)  
       % +1 edge dof
       NE=NE+1;
       % update edge_connectivity
       edge_connectivity( n1, n2 )=NE;
       edge_connectivity( n1, n2 )=NE;
       % update edge with the node dofs 
       edge(NE,[1,2])=edges(n,[1,2]);
   end
   
   NE_tmp=edge_connectivity( n1, n2 );
      
   if(isempty(E_to_T{NE_tmp}))
   E_to_T{NE_tmp}=cell(1,1);
   E_to_T{NE_tmp}{end}=T;
   else
   E_to_T{NE_tmp}{end+1}=T;   
   end
   % create a structure like elem, but with edge dofs instead of vertex dofs 
   elemE(T,n)=NE_tmp; 
   end

   for n=1:node_per_elem
   % node dof of the element
   vv=elem(T,n);
   if(isempty(N_to_T{vv}))
   N_to_T{vv}=cell(1,1);
   N_to_T{vv}{end}=T;
   else
   N_to_T{vv}{end+1}=T;   
   end       
   end
   
end

 
 E_to_T=E_to_T(~cellfun('isempty',E_to_T))  ;
end

