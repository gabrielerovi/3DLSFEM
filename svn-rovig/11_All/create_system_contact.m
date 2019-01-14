
function [Ant,Complementarity,bnt,bnt1,Pnt]=create_system_contact(parameters,mesh,h,Patch_Node_Monotone)
qrule=parameters.qrule;
L=length(mesh);
toll=10^-12;

P.P1CtoP1F=P1CtoP1F3D(mesh,Patch_Node_Monotone);
P.RTCtoRTF= RTCtoRTF3D(mesh,qrule);

if(L>1)
for lev=1:L-1
[m,n]=size(P.RTCtoRTF{lev});
P.RTCtoRTFss{lev}= [P.RTCtoRTF{lev}   sparse(m,n)        sparse(m,n);
                    sparse(m,n)       P.RTCtoRTF{lev}  sparse(m,n);
                    sparse(m,n)       sparse(m,n)       P.RTCtoRTF{lev}];
[m,n]=size(P.P1CtoP1F{lev});
P.P1CtoP1Fuu{lev}= [P.P1CtoP1F{lev}, sparse(m,n),       sparse(m,n);
                    sparse(m,n),     P.P1CtoP1F{lev},   sparse(m,n);
                    sparse(m,n),     sparse(m,n),       P.P1CtoP1F{lev}];
end
end

Pnt=[];
if(parameters.L>1)
 Pnt=ArnoldPnt(P,mesh);
end


[Ant,Complementarity,bnt,bnt1]=create_system_LSelasticityContact(parameters,mesh,h,P);

end



