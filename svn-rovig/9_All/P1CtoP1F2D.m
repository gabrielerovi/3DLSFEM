function [P] =P1CtoP1F2D(mesh,qrule)
% Compute the L2 projection between P1C and P1F
% M_FF x = M_FC y
% x = M_FF^(-1) M_FC y
% x = T y, with T=L2 projection

L=size(mesh);
L=L(1);
node_per_elem=mesh{1}.node_per_elem;
NDCtoNDF=cell(L-1,1);
massNDFNDC=cell(L-1,1);
massNDF=cell(L-1,1);
P=cell(L-1,1);
q_point=zeros(1,1);




%loop on each connection between subsequent levels 
for lev=2:L
F=lev;
C=lev-1;
NF=mesh{F}.N;
NC=mesh{C}.N;
NTF=mesh{F}.NT;
NTC=mesh{C}.NT;

% massNDFNDC: NDC -> NDF
%massNDFNDC{C}=zeros(NF,NC);
massNDFNDC{C}=sparse(NF,NC);
% massNDFNDC: NDF -> NDF
%massNDF{C}=zeros(NF,NF);
massNDF{C}=sparse(NF,NF);
%loop on the elements of the fine grid
for tF=1:NTF
    % nodes of the fine element
    elemF=mesh{F}.elem(tF,:);
    nodeF=mesh{F}.node(elemF,:);
    [q_point,weights,area]=quadrature_points_2D(qrule,nodeF);
    number_of_qp=length(q_point(:,1));
    P1_basisF = phiP12D(q_point,nodeF);
    
    %%%%%%% compute M_FF  %%%%%%%
    for nF1=1:node_per_elem
        
        for nF2=1:node_per_elem
                 value=zeros(number_of_qp,1);
                 for qp=1:number_of_qp
                     value(qp)=weights(qp)*P1_basisF(nF1,qp) * P1_basisF(nF2,qp);
                 end
                 q=sum(value)*area;
                 raw=elemF(nF1);
                 col=elemF(nF2);
                 massNDF{C}(raw,col)=massNDF{C}(raw,col)+q;            
        end
    end
    
    %%%%%%% compute M_FC %%%%%%%
    % from the fine element, recover the coarse element that contains this edge
    tC=ceil(tF/4); 
    elemC=mesh{C}.elem(tC,:);
    nodeC=mesh{C}.node(elemC,:);
    P1_basisC = phiP12D(q_point,nodeC);
    for nF=1:node_per_elem        
        for nC=1:node_per_elem
                 value=zeros(number_of_qp,1);
                 for qp=1:number_of_qp
                     value(qp)=weights(qp)*P1_basisF(nF,qp) * P1_basisC(nC,qp);
                 end
                 q=sum(value)*area;
                 raw=elemF(nF);
                 col=elemC(nC);
                 massNDFNDC{C}(raw,col)=massNDFNDC{C}(raw,col)+q;            
        end
    end
end
P{C}= (massNDF{C})^(-1) * massNDFNDC{C};
end
end