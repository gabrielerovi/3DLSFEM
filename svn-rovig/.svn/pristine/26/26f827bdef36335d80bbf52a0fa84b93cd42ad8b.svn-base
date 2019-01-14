function P1toRT=P1_to_RT(mesh)

L=size(mesh);
L=L(1);
P1toRT=cell(L,1);

for lev=1:L
N=mesh{lev}.N;     
NE=mesh{lev}.NE;    
%P1toRT{lev}=zeros(NE,N);
P1toRT{lev}=sparse(NE,N);
edge=mesh{lev}.edge;
for ee=1:NE
    side=sort(edge(ee,:));
    P1toRT{lev}(ee,side(1))=-1;
    P1toRT{lev}(ee,side(2))=+1;
    
end

end

end