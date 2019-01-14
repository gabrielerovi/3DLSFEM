function [x_k,activesetold]=activesetResidual(A,b,c,activeset)
 
toll=10^(-64);

% A system matrix
% b system rhs
% x system solution
% c vector containing the constraints for the whole x
% in case x_i is unconstrained, c_i=Inf
% w workingset

% we solve for the system A x + lambda =f
% subject to the constraint x<=c

res_k=1;

while(1)
A_k=A;
b_k=b;

% fix the active set
for jj=1:length(activeset)
ii=activeset(jj);    
A_k(ii,:)=0;
A_k(ii,ii)= 1;
b_k(ii)= c(ii);
end

x_k=A_k\b_k;
res_k=b_k-A_k*x_k;

activesetold=activeset;

inactiveset=find(b_k-A*x_k<toll)';
activeset=unique(setdiff(activeset,inactiveset));
tmp=find(x_k-c>=toll)';
activeset=unique([activeset,tmp]);
% inactiveset=find(b_k-A_k*x_k<0);
% activeset=unique(setdiff(activeset,inactiveset));

if(isequal(activeset,activesetold)) 
    break
end
end


end