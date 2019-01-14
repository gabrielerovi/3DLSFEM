function [x,cont,normresvec]= LS_PC(A11,A12,A22,b1,b2,P1,P2,toll,MaxIT)


A= [ A11  A12 ; 
     A12' A22 ];
b=[b1;b2]; 

n1= length(A11(:,1)); 
n2= length(A22(:,1));

first=1:1:n1;
second=n1+1:1:n1+n2;

x=zeros(n1+n2,1);

res=[b1;b2];

z1=P1\res(first);
z2=P2\res(second);
z=[z1;z2];
p=z;
cont=0;

while( cont<MaxIT)
    cont=cont+1;
    Ap=A*p;
    pAp=p'*Ap;
    zres=z'*res;
    
    alpha = zres / pAp;
    
    x= x + alpha * p;
    
    res=res-alpha*Ap;
    normres=norm(res);
    normresvec(cont)=normres;
    if(normres<toll)
        break
    end
    znew1= P1\ res(first);
    znew2= P2\ res(second);
    z=[znew1;znew2];
    
    beta= z'*res/ (zres);
    p=z+beta*p;
end



end