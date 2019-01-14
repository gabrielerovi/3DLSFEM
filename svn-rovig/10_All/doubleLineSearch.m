function [d]=doubleLineSearch (A11,A12,A21,A22,correction1, correction2, res1,res2)

toll=10^(-10);
b(1,1)=correction1'*res1;
b(2,1)=correction2'*res2;

M(1,1)= correction1'*A11*correction1;
M(1,2)= 0.5 * (correction1'*A12*correction2 + correction2'*A21*correction1);
M(2,1)= M(1,2);
M(2,2)= correction2'*A22*correction2;

if(det(M)<toll)
d=[1;1];
else
d=M\b;
end


end