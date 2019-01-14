function P=P1CtoP1F(mesh,Patch_Node_Monotone)

L=length(mesh);
P=[];
for lev=1:L-1
    
    PNM=Patch_Node_Monotone{lev};
    NC=mesh{lev}.N;
    NF=mesh{lev+1}.N;
    P{lev}=sparse(NF,NC);

for nC1=1:NC
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Number of nodes belonging to the coarse patch of nC  %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for neighb=1:length(PNM{nC1})
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%  nC1 ----- nF ----- nC2  %%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        nF=PNM{nC1}{neighb}(2);
        nC2=PNM{nC1}{neighb}(1);
        P{lev}(nC1,nC1)=1;
        P{lev}(nC2,nC2)=1;
        P{lev}(nF,nC1)=0.5;
        P{lev}(nF,nC2)=0.5;
    end    
  end
end


end