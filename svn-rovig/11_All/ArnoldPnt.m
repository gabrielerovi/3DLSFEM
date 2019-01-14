function Pnt=ArnoldPnt(P,mesh)

L=length(mesh);

for lev=1:L

[M_Normal_Tangent{lev}] = MatrixOnGammaCwithNormalTangentComponents(mesh{lev});

end


dim=length(mesh{1}.node(1,:));

for lev=1:L-1
    
    F=lev+1;
    C=lev;
    
    NFF=mesh{F}.NF;
    NF=mesh{F}.N;
    NFC=mesh{C}.NF;
    NC=mesh{C}.N;
    
    Ptot{C}=sparse(dim*(NFF+NF),dim*(NFF+NF));
    
    Ptot{C}= [ P.RTCtoRTFss{C}           sparse(dim*NFF,dim*NC );
               sparse(dim*NF,dim*NFC )   P.P1CtoP1Fuu{C}];
    
     Pnt{C} = M_Normal_Tangent{F} * Ptot{C} * M_Normal_Tangent{C}';
end

end
