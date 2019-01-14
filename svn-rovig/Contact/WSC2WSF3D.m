function WorkingSetF= WSC2WSF3D(WorkingSetC,mesh_parameters)

NFC=mesh_parameters.NFC;
NFF=mesh_parameters.NFF;
NC=mesh_parameters.NC;
NF=mesh_parameters.NF;

% We transfer the WSC to WSF

if(length(WorkingSetC)<3*NFC+3*NC)
    tmp=sparse(3*NFC+3*NC,1);
    tmp(WorkingSetC)=1;
    WorkingSetC=tmp;
end

WorkingSetF=sparse(3*NFF+3*NF,1);

tmp= find(WorkingSetC(3*NFC+1:3*NFC+NC)==1);
WorkingSetF( 3* NFF  + tmp)=1;

facesC=find(WorkingSetC(1:NFC)==1);
for ii=1:length(facesC)
    fC=facesC(ii);
    faceF=mesh_parameters.Patch_Face_Monotone{fC};    
    WorkingSetF(faceF)=1;
end


end