
function [Ant,Complementarity,bnt,bnt1,Pnt]=create_system_contact(parameters,mesh,h,Patch_Node_Monotone)
qrule=parameters.qrule;
L=length(mesh);
toll=10^-12;
% create projections
P.RTCtoRTF= RTCtoRTF2D(mesh,qrule);
%P.P1CtoP1F =P1CtoP1F2D(mesh,qrule);
P.P1CtoP1F=P1CtoP1F(mesh,Patch_Node_Monotone);
Antbc=[];
if(L>1)
for lev=1:L-1
% here we remove errors in the computation of P1 projection    
P.P1CtoP1F{lev}(abs(P.P1CtoP1F{lev})<toll)=0;

[m,n]=size(P.RTCtoRTF{lev});
P.RTCtoRTFss{lev}= [P.RTCtoRTF{lev} zeros(m,n);
                    zeros(m,n)       P.RTCtoRTF{lev}];
[m,n]=size(P.P1CtoP1F{lev});
P.P1CtoP1Fuu{lev}= [P.P1CtoP1F{lev} zeros(m,n);
                    zeros(m,n)      P.P1CtoP1F{lev}];
end
end




 % make the mesh circular
if(parameters.create_only_square==0)
%     mesh=mesh_make_it_circle(mesh);
end


Pnt=[];
if(parameters.L>1)
 Pnt=ArnoldPnt(P,mesh);
end




 
[Ant,Complementarity,bnt,bnt1]=create_system_LSelasticityContact(parameters,mesh,h,P);

end



