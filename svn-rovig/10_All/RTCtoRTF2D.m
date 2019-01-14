function RTCtoRTF= RTCtoRTF2D(mesh,qrule)

edge_per_elem=mesh{1}.edge_per_elem;
epsilon=0.0000001;
L=size(mesh);
L=L(1);
if(L<2)
 NDtoRT=[];
 RTCtoRTF=[];
else
NDtoRT=cell(L-1,1);
RTCtoRTF=cell(L-1,1);
% kk gives the opposite vertex of a given edge
kk(1)=3; kk(2)=1; kk(3)=2; 
%  3\
%  | \
%  c  b
%  |   \
%  1--a--2
%  
%  edge a=1,edge b=2, edge c=3
%





% we loop on all intergrid connections (L-1)
for lev=2:L
    
    F=lev;
    C=lev-1;
    NF=mesh{lev}.N;
    NEF=mesh{F}.NE;
    NC=mesh{C}.N;
    NEC=mesh{C}.NE;
    
    %RTCtoRTF{C}=zeros(NEF,NEC);
    RTCtoRTF{C}=sparse(NEF,NEC);
    
% we loop on each fine edge, that we need only for the domain and the
% normal of the integral
for eF=1:NEF 
    edge=mesh{F}.node(mesh{F}.edge(eF,:), : );
    normal=[edge(2,2)-edge(1,2),   -edge(2,1)+edge(1,1)];
    normal=normal/norm(normal);
    % fine element that share this edge
    TF=mesh{F}.E_to_T{eF};
    NT=size(TF);NT=NT(2);  
    for t=1:NT
    % from the fine element, recover the coarse element that contains this edge
    TC=ceil(TF{t}/4); 
    % define the nodes of the coarse element
    nodeC=mesh{C}.node(mesh{C}.elem(TC,:),:);
    % loop on the edges of the coarse element TC
    for eC=1:edge_per_elem
        pC=nodeC(kk(eC),:);       
        raw=eF;
        col=mesh{C}.elemE(TC,eC);
       if(abs(RTCtoRTF{C}(raw,col))<epsilon)
       RTCtoRTF{C}(raw,col)=edge_integral_RTF_RTC2D(eC,edge,nodeC,qrule,normal,pC);
       end
    end
   end
end
end
end

end