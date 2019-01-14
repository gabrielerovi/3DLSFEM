function WorkingSetF= WSC2WSF(WorkingSetC,mesh,maps,C,F)

NEC=mesh{C}.NE;
NEF=mesh{F}.NE;
NC=mesh{C}.N;
NF=mesh{F}.N;

% We transfer the WSC to WSF

if(length(WorkingSetC)<2*NEC+2*NC)
    tmp=sparse(2*NEC+2*NC,1);
    tmp(WorkingSetC)=1;
    WorkingSetC=tmp;
end

WorkingSetF=sparse(2*NEF+2*NF,1);

tmp= find(WorkingSetC(2*NEC+1:2*NEC+NC)==1);
WorkingSetF( 2 * NEF  + tmp)=1;

edgesC=find(WorkingSetC(1:NEC)==1);
for ii=1:length(edgesC)
    eC=edgesC(ii);
    edgeF=maps.Patch_Edge_Monotone{C}{eC};    
    WorkingSetF(edgeF)=1;
end


end