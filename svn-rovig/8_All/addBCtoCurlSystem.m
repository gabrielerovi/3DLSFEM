function CC=addBCtoCurlSystem(C,NeumannNode,BCneighb,BCneighbCoeff)

M=length(NeumannNode);
N=length(C(1,:))/2;
CC=C;
for nn=1:M
    
    ii=NeumannNode(nn);
    jj=ii+N;
    BCii=BCneighb(nn,:);
    BCjj=BCii+N;
    BCcoeff=BCneighbCoeff{nn};
    CC(ii,:)=0;
    CC(jj,:)=0;
    CC(ii,ii)=BCcoeff(1);
    CC(jj,jj)=BCcoeff(1);
    
    CC(ii,BCii(1))=BCcoeff(2);
    CC(ii,BCii(2))=BCcoeff(3);
    
    CC(jj,BCjj(1))=BCcoeff(2);
    CC(jj,BCjj(2))=BCcoeff(3);
    
    
end



end