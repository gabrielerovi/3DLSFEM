function WorkingSetF= WSC2WSF(WorkingSetC,mesh_parameters)

NEC=mesh_parameters.NEC;
NEF=mesh_parameters.NEF;
NC=mesh_parameters.NC;
NF=mesh_parameters.NF;

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
    edgeF=mesh_parameters.Patch_Edge_Monotone{eC};    
    WorkingSetF(edgeF)=1;
end


end