function [xn,H]=HouseHolderTransformation(x,normal)
toll=10^(-10);
dim=length(x);
% normalize normal, just to be sure
normal=normal/norm(normal);
% define first basis vector
E1=zeros(dim,1);
E1(1)=1;
% define v=0.5*(normal-E1)
v=0.5*(normal-E1);
% normalize vector
if(norm(v)>toll)
v=v/norm(v);
end
sizev=size(v);
% HouseHolderTransformation
if(sizev(1)==1)
    v=v';
end
H=eye(dim)-2*v*v';

xn=H*x;



end