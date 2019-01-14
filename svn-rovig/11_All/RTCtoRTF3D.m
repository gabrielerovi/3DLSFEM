function RTCtoRTF= RTCtoRTF3D(mesh,qrule)

face_per_elem=mesh{1}.face_per_elem;
epsilon=0.0000001;
L=size(mesh);
L=L(1);
if(L<2)
 RTCtoRTF=[];
else
RTCtoRTF=cell(L-1,1);
% kk gives the opposite vertex of a given face
kk(1)=1; kk(2)=2; kk(3)=3; kk(4)=4;


% we loop on all intergrid conNFCtions (L-1)
for lev=2:L
    
    F=lev;
    C=lev-1;
    NF=mesh{lev}.N;
    NFF=mesh{F}.NF;
    NC=mesh{C}.N;
    NFC=mesh{C}.NF;
    
    RTCtoRTF{C}=sparse(NFF,NFC);
    
% we loop on each fine face, that we need only for the domain and the
% normal of the integral
for fF=1:NFF 
    nodes_coord_of_faceF=mesh{F}.node(mesh{F}.face(fF,:), : );
    normalF =RT03Dnormal(nodes_coord_of_faceF);
    % fine elemFnt that share this face
    TF=mesh{F}.F_to_T{fF};
    NT=size(TF);NT=NT(2);  
    for t=1:NT
    % from the fine elemFnt, recover the coarse element that contains this face
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% QUI CONTROLLA CHE SIA EFFETTIVAMENTE COSI CAZZO %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    TC=ceil(TF{t}/8); 
    % dfFine the nodes of the coarse elemFnt
    nodes_coord_of_tetrahedronC=mesh{C}.node(mesh{C}.elem(TC,:),:);
    % loop on the faces of the coarse elemFnt TC
    for fC=1:face_per_elem
        pC=nodes_coord_of_tetrahedronC(kk(fC),:);       
        row=fF;
        col=mesh{C}.elemF(TC,fC);
       if(abs(RTCtoRTF{C}(row,col))<epsilon)
       RTCtoRTF{C}(row,col)=face_integral_RTF_RTC3D(fC,nodes_coord_of_faceF,nodes_coord_of_tetrahedronC,qrule,normalF,pC);
       end
    end
   end
end
end
end

end