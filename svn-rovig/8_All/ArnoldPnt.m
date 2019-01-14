function Pnt=ArnoldPnt(P,mesh)

L=length(mesh);
householder=1;

for lev=1:L

[M_Normal_Tangent{lev},M_Normal_TangentT{lev}] = MatrixOnGammaCwithNormalTangentComponents(mesh{lev},householder);

end


for lev=1:L-1
    
    F=lev+1;
    C=lev;
    
    NEF=mesh{F}.NE;
    NF=mesh{F}.N;
    NEC=mesh{C}.NE;
    NC=mesh{C}.N;
    
    Ptot{C}=sparse(2*NEF+2*NF,2*NEC+2*NC );
    
    Ptot{C}= [ P.RTCtoRTFss{C}    sparse(2*NEF,2*NC );
                 sparse(2*NF,2*NEC )   P.P1CtoP1Fuu{C}];
    
     Pnt{C} = M_Normal_Tangent{F} * Ptot{C} * M_Normal_Tangent{C}';
%     Pnt.RTCtoRTF{lev}=  M_Normal_Tangent{lev+1} * P.RTCtoRTF{lev} * M_Normal_Tangent{lev}';
%     Pnt.P1CtoP1F{lev}=  M_Normal_Tangent{lev+1} * P.P1CtoP1F{lev} * M_Normal_Tangent{lev}';
end

end
