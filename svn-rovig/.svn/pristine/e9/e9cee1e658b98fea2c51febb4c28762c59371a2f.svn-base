    function [C11_lev,C12_lev,C22_lev,NDtoRT]=assemblingC(C11_lev,C12_lev,C22_lev,NDtoRT,mesh)


for lev=1:length(mesh)
N=mesh{lev}.N;
reorder=mesh{lev}.ND_nodes_reorder;
Lengths=mesh{lev}.ND_lengths;
edge_boundary=mesh{lev}.ND_nodes_on_edge_boundary;

C11_lev{lev}=C11_lev{lev}(reorder,reorder);
C12_lev{lev}=C12_lev{lev}(reorder,reorder);
C22_lev{lev}=C22_lev{lev}(reorder,reorder);
NDtoRT{lev}=NDtoRT{lev}(:,reorder);

% loop on the first boundary on which we impose potential=0 
for nn=1:Lengths(1)
    C11_lev{lev}(nn,:)=zeros(1,N);
    C11_lev{lev}(nn,nn)=1;   
    C12_lev{lev}(nn,:)=zeros(1,N);
    C12_lev{lev}(nn,nn)=1;   
    C22_lev{lev}(nn,:)=zeros(1,N);
    C22_lev{lev}(nn,nn)=1;   
end
% loop on the other boundary on which we impose potential=cost
% so the first node is internal and there the potential is computed
% normally. once the value is known, it is transferred to the other nodes
% of that boundary. so, for a boundary, we impose dirichlet on all the
% nodes, with the exception of the first one
for nn1=2: (length(Lengths)-1)
    
    S=sum(Lengths(1:nn1-1));
    
    for nn2=2:Lengths(nn1)
    C11_lev{lev}(S+nn2,:)=zeros(1,N);
    C11_lev{lev}(S+nn2,S+nn2)=1; 
    C12_lev{lev}(S+nn2,:)=zeros(1,N);
    C12_lev{lev}(S+nn2,S+nn2)=1; 
    C22_lev{lev}(S+nn2,:)=zeros(1,N);
    C22_lev{lev}(S+nn2,S+nn2)=1;     
    end
    
end





end


end
